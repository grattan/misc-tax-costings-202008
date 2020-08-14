
library(expm) # %^%
library(data.table)
library(taxstats)
library(grattan)
library(hutilscpp)
library(magrittr)
library(hutils)

cgt_reducer_due_to_policy <- 0.015
new.discount <- 0.25



s1718 <-
  fread("../taxstats1718/2018_sample_file.csv") %>%
  .[taxstats::age_range_decoder, age_range_description := i.age_range_description, on = "age_range"] %>%
  .[, min_age := if_else(grepl("to", age_range_description),
                         as.numeric(stringr::str_extract(age_range_description, "^[0-9]{2}")),
                         if_else(grepl("70", age_range_description),
                               70,
                               15)),
  keyby = "age_range"] %>%
  .[, max_age := min_age + 5] %>%
  .[, age_imp := runif(.N, min_age, max_age), keyby = "age_range"] %>%
  .[, Tot_inc_amt := Tot_IncLoss_amt] %>%
  .[, WEIGHT := 50] %>%
  .[]



revenue_quarantine_NG_2 <- function(.target_sample_file,
                                    fy.year,
                                    new_CGT_discount = new.discount,
                                    age_based = TRUE,
                                    Z = 20) {

  cgt_discount_distribution <-
    .target_sample_file %>%
    .[Tot_CY_CG_amt > 0] %>%
    .[, apparent_discount := round((1 - Net_CG_amt / Tot_CY_CG_amt) * 4) / 4] %>%
    .[, .(value = sum(Tot_CY_CG_amt)), keyby = .(apparent_discount)] %>%

    .[, density := value / sum(value)] %>%
    .[order(apparent_discount)] %>%
    .[, cumdensity := cumsum(density)]

  probNoDiscount <-
    cgt_discount_distribution %>%
    .[apparent_discount == 0] %>%
    .[["cumdensity"]]

  target_sample_file <- copy(.target_sample_file)

  probCG_by_age <-
    target_sample_file %>%
    .[Gross_rent_amt > 0] %>%
    .[, .(probCG_overall = mean(Net_CG_amt > 0),
          # after discount
          average_CG = 1.5 * mean(Net_CG_amt)),
      keyby = .(age = floor(age_imp))]

  cg_by_age <-
    target_sample_file %>%
    .[Net_CG_amt > 0 | Gross_rent_amt > 0] %>%
    .[, cgt_up := (1 - cgt_reducer_due_to_policy)] %>%
    .[, Tot_CG_after_policy := cgt_up * Tot_CY_CG_amt] %>%
    .[, has_CG := Net_CG_amt > 0] %>%
    .[, .(probCG_overall = mean(has_CG),
          total_cg = sum(Tot_CG_after_policy),
          expected_cg = mean(Tot_CG_after_policy),
          expected_tot_cg_given_cg_occurs = mean(Tot_CG_after_policy[Net_CG_amt > 0])),
      keyby = (age = floor(age_imp))]

  if (!age_based) {
    # OK if constant


    M4 <- matrix(0, nrow = 2*Z, ncol = 2*Z)
    M4[row(M4) == col(M4) - 2] <- probNG_if_investor
    M4[row(M4) == col(M4) - 1] <- 1 - probNG_if_investor
    M4[row(M4) %% 2 == 0] <- 0
    M4[row(M4) %% 2 == 0 & row(M4) == col(M4)] <- 1
    M4[2*Z - 1, 2*Z] <- 1

    stopifnot(all(rowSums(M4) == 1))

    P <- M4 %^% (Z + 1)

    init <- c(1, rep(0, Z - 1))

    # even entries of P are the absorbing states
    E <- init %*% P[c(TRUE, FALSE),c(FALSE,TRUE)]

    accumulator <-
      data.table(rand = cumsum(E[1,]),
                 loss_multiplier = 1:(Z)) %>%
      setkey(rand)

    forecast_loss_multiple <- function(.data) {
      setkey(.data, rand) %>%
        accumulator[., roll = -Inf]
    }

  } else {

    # age based-method
    # q ~ probNG_if_investor
    q_by_age <-
      copy(s1718) %>%
      .[Gross_rent_amt > 0] %>%
      .[, .(q = mean(Net_rent_amt < 0) * (1 - mean(Net_CG_amt > 0))),
        keyby = "age_range"]

    q_by_integer_age <-
      data.table(age_range = rep(q_by_age$age_range, each = 5)) %>%
      merge(q_by_age, by = "age_range") %>%
      merge(age_range_decoder) %>%
      .[, min_age := if_else(grepl("to", age_range_description),
                             as.numeric(stringr::str_extract(age_range_description, "^[0-9]{2}")),
                             if_else(grepl("70", age_range_description),
                                     70,
                                     15))] %>%
      .[, age := min_age - 1 + seq_len(.N), keyby = "age_range"] %>%
      setkey(age) %>%
      .[] %>%
      .[, .(age, q)]

    # We want a wide data frame which shows, for each
    # age, the probability of going from n_ahead to n_ahead + 1
    # so [i, j] represents the probability that a person aged i
    # who *will* negative gear for at least j years will negatively gear
    # in the year j + 1.

    # All the information to produce this dataframe/matrix is contained in q
    # so we just need to extract it.
    # We keep applying lead until just before there will be a column of pure NA.
    j <- 1
    setnames(q_by_integer_age, "q", "q1")
    alloc.col(q_by_integer_age)
    while (TRUE) {
      j <- j + 1
      set(q_by_integer_age,
          j = paste0("q", j),
          value = shift(.subset2(q_by_integer_age, paste0("q", j - 1L)),
                        type = "lead"))

      # breaker:
      if (any(colSums(is.na(q_by_integer_age)) == nrow(q_by_integer_age) - 1)) {
        # We're done. Before we `break;` let's clean up.
        # We assume that people persist in this model beyond 75
        na_locf <- function(x) {
          zoo::na.locf(x, na.rm = FALSE)
        }

        for (j in seq_along(q_by_integer_age)) {
          set(q_by_integer_age, j = j, value = na_locf(.subset2(q_by_integer_age, j)))
        }
        break
      }
    }
    setDT(q_by_integer_age)

    q_by_integer_age_long <-
      q_by_integer_age %>%
      melt.data.table(id.vars = c("age"),
                      variable.name = "n_ahead",
                      value.name = "prob_NG_if_investor") %>%
      .[, n_ahead := as.numeric(gsub("q", "", n_ahead))] %>%
      setkey(age, n_ahead) %>%
      .[]

    E_by_age <- function(the_age) {
      probNG_if_investor <-
        q_by_integer_age_long %>%
        .[age == the_age] %$%
        prob_NG_if_investor

      Z <- length(probNG_if_investor)

      # Construct a matrix as before:
      # odd states represent continuing negative gearing
      # even states are absorbing
      #
      # The idea is that once an individual stop negative gearing
      # that person never negative gears and their accumulated losses
      # are those thitherto accumulated.
      MA <- matrix(0, nrow = 2*Z, ncol = 2*Z)
      # if still NG, advance
      MA[row(MA) == col(MA) - 2] <- rep(probNG_if_investor, each = 2)[-c(1, 2*Z)]  ## dropped indices fall off the 'north' and 'east' the matrix
      # else fall into absorbing state
      MA[row(MA) == col(MA) - 1] <- rep(1 - probNG_if_investor, each = 2)[-c(2*Z)]

      # even entires are absorbing
      MA[row(MA) %% 2 == 0] <- 0
      MA[row(MA) %% 2 == 0 & row(MA) == col(MA)] <- 1

      # Final state hack to ensure stochastic.
      # We take the last state as absorbing, even if the person is (by some miracle) still NG.
      MA[2*Z - 1, 2*Z] <- 1

      # MA must be stochastic
      stopifnot(all(rowSums(MA) == 1))

      P <- MA %^% (Z + 1)

      init <- c(1, rep(0, Z - 1))

      # even entries of P are the absorbing states
      E <- init %*% P[c(TRUE, FALSE), c(FALSE,TRUE)]

      data.table(
        age = rep(the_age, Z),
        loss_multiplier = 1:Z,
        E = as.numeric(t(E))
      )
    }

    accumulator_by_age <-
      lapply(unique(q_by_integer_age_long$age), E_by_age) %>%
      rbindlist(.) %>%
      .[, rand := cumsum(E), keyby = "age"] %>%
      setkey(age, rand)

    forecast_loss_multiple <- function(.data) {
      setkey(.data, age, rand) %>%
        accumulator_by_age[., roll = -Inf]
    }

    probCG__by_age <-
      target_sample_file %>%
      .[, .(prob__CG_event = mean(Tot_CY_CG_amt > 0)),
        keyby = (age_at_Z = floor(age_imp))] %>%
      setkey(age_at_Z)

    forecast_cg_event <- function(.data){
      setkey(.data, age_at_Z) %>%
        merge(probCG__by_age)
    }

  }

  target_sample_file %>%
    selector(Net_CG_amt, Sw_amt, ETP_txbl_amt,
             Tot_ded_amt, NPP_loss_claimed, PP_loss_claimed,
             Alow_ben_amt, Tot_inc_amt, Net_rent_amt, Taxable_Income, age_imp,
             WEIGHT) %>%
    # mutate_each(funs(round_cols)) %>%

    # We need to work out how much the discount will change Net_CG_amt.
    # We use the density of the apparent_discount to randomly
    # assign a new discount rate.  Note that this does not
    # adequately take into account the impact of capital losses

    .[, new_Net_CG_amt := (1 - cgt_reducer_due_to_policy) * Net_CG_amt * if_else(runif(.N) < probNoDiscount,
                                                                                 1,
                                                                                 (1 - new_CGT_discount) / 0.50)] %>%

    # Excluding CG for the moment, we create a new total income amount
    # which quarantines salary from NG.  Rental losses may be deducted
    # against investment income (including capital gains).  However,
    # we assume that they are deducted against capital gains only
    # in the future

    # Other_Income means income (possibly -ve) that is
    # not salary, CG, or rental
    .[, quarantined_income := Sw_amt + Alow_ben_amt + ETP_txbl_amt] %>%
    .[, Other_Income := Tot_inc_amt - quarantined_income - Net_CG_amt - Net_rent_amt] %>%
    .[, new_Tot_inc_amt_no_CG := quarantined_income + if_else(Other_Income <= 0,
                                                              # Nothing for rent to deduct against,
                                                              # but Other_Income is permitted to be -ve
                                                              Other_Income,
                                                              # Net_rent is negative and is allowed to
                                                              # reduce Other income till exhausted
                                                              pmax0(Other_Income + Net_rent_amt))] %>%
    .[,
      # should be positive
      accumulated_losses_year_0 := -1 * if_else(Net_rent_amt >= 0,
                                                0,
                                                as.double(pmin0(pmax0(Other_Income) + Net_rent_amt)))
      # Not required.
      # ,new_Tot_inc_amt = new_Tot_inc_amt_no_CG + Net_CG_amt
      ] %>%
    .[, rand := runif(.N)] %>%
    .[, age := floor(age_imp)] %>%
    forecast_loss_multiple %>%


    .[, accumulated_losses_year_Z := accumulated_losses_year_0 * loss_multiplier] %>%
    .[, age_at_Z := age + loss_multiplier] %>%
    merge(cg_by_age, by.x = "age_at_Z", by.y = "age") %>%
    forecast_cg_event %>%


    .[, income_reduction := if_else(rand < prob__CG_event & accumulated_losses_year_Z > 0,
                                    # in this case we expect a capital gain
                                    # and the concommitant income reduction is
                                    # ither the gain (if the accumulated losses exceed it)
                                    # or the losses themselves
                                    pmaxV(expected_tot_cg_given_cg_occurs, accumulated_losses_year_Z),
                                    accumulated_losses_year_Z)] %>%
    .[, income_reduction := pminC(accumulated_losses_year_Z, 40e3)] %>%
    .[, gain_residual := if_else(rand < prob__CG_event & accumulated_losses_year_Z > 0,
                                 pmaxC(expected_tot_cg_given_cg_occurs - accumulated_losses_year_Z, 0),
                                 0)] %>%

    # This is the standard calculation for taxable income
    .[, new_Taxable_Income_yearZ := pmax0(pmax0(Other_Income + new_Net_CG_amt + ((1 - new_CGT_discount) / 0.50) * gain_residual - income_reduction) + quarantined_income - Tot_ded_amt - NPP_loss_claimed - PP_loss_claimed)] %>%
    .[, new_Taxable_Income_year0 := pmax0(pmax0(Other_Income + new_Net_CG_amt) + quarantined_income - Tot_ded_amt - NPP_loss_claimed - PP_loss_claimed)] %>%

    # pmaxC redundant because taxstats have already done it
    .[, prev_tax := income_tax(pmaxC(Taxable_Income, 0), fy.year = fy.year)] %>%
    .[, post_tax_yearZ := income_tax(pmaxC(new_Taxable_Income_yearZ, 0), fy.year = fy.year)] %>%
    .[, post_tax_year0 := income_tax(new_Taxable_Income_year0, fy.year = fy.year)] %>%
    .[, minimum_taxable_income := pmax0(quarantined_income + new_Net_CG_amt - Tot_ded_amt - NPP_loss_claimed - PP_loss_claimed)] %>%
    .[, tax_increase := post_tax_yearZ - prev_tax] %>%
    .[, tax_increase_year0 := post_tax_year0 - prev_tax] %>%
    .[, residual_gain_effect := income_tax(pmaxC(pmaxC(Other_Income + new_Net_CG_amt - income_reduction,
                                                       0) + quarantined_income - Tot_ded_amt - NPP_loss_claimed - PP_loss_claimed,
                                                 0),
                                           fy.year = fy.year) - post_tax_yearZ] %>%
    .[]
}


revenue_CGT_discount <- function(.target_sample_file = target_sample_file,
                                 fy.year = "2019-20",
                                 CGT.discount = new.discount,
                                 .old.discount = 0.50) {
  cgt_discount_distribution <-
    .target_sample_file %>%
    .[Tot_CY_CG_amt > 0] %>%
    .[, apparent_discount := round((1 - Net_CG_amt / Tot_CY_CG_amt) * 4) / 4] %>%
    .[, .(value = sum(Tot_CY_CG_amt)), keyby = .(apparent_discount)] %>%

    .[, density := value / sum(value)] %>%
    .[order(apparent_discount)] %>%
    .[, cumdensity := cumsum(density)]

  probNoDiscount <-
    cgt_discount_distribution %>%
    .[apparent_discount == 0] %>%
    .[["cumdensity"]]

  .target_sample_file %>%
    copy %>%
    # Assume 2 % fall in gains.
    .[, new_Net_CG_amt := (1 - cgt_reducer_due_to_policy) * if_else(runif(.N) < probNoDiscount,
                                                                    as.double(Net_CG_amt),
                                                                    (Net_CG_amt / (1 - .old.discount)) * (1 - CGT.discount))] %>%
    .[, new_Tot_inc_amt := Tot_inc_amt - Net_CG_amt + new_Net_CG_amt] %>%
    .[, new_Tot_ded_amt := Tot_ded_amt + 0] %>%   # no change
    .[, new_Taxable_Income := pmax0(Taxable_Income - Net_CG_amt + new_Net_CG_amt)] %>%
    .[, old_tax := income_tax(Taxable_Income, fy.year = fy.year, age = age_imp)] %>%
    .[, new_tax := income_tax(new_Taxable_Income, fy.year = fy.year, age = age_imp)] %>%
    .[, tax_increase := new_tax - old_tax]
}



revenue_quarantine_NG_2(copy(s1718), "2019-20", new_CGT_discount = 0.25) %>%
  .[, grattan:::print.revenue_foregone(sum(tax_increase * WEIGHT))]

revenue_CGT_discount(copy(s1718), "2019-20",CGT.discount = 0.25) %>%
  .[, grattan:::print.revenue_foregone(sum(tax_increase * WEIGHT))]
