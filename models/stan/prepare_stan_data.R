#' Prepare LL/SS trial arrays from sujeto_*.csv (same rules as data/parse_1.m
#' + models/prepareIntertemporalChoiceData.m).

prepare_stan_data <- function(data_dir, mu_prec_scale = 1) {
  n_participants <- 25L
  n_trials <- 240L
  reward_a <- matrix(NA_real_, n_participants, n_trials)
  reward_b <- matrix(NA_real_, n_participants, n_trials)
  time_a <- matrix(NA_real_, n_participants, n_trials)
  time_b <- matrix(NA_real_, n_participants, n_trials)
  choice <- matrix(NA_integer_, n_participants, n_trials)

  for (i in seq_len(n_participants)) {
    f <- file.path(data_dir, sprintf("sujeto_%d_tiempo.csv", i))
    T <- utils::read.csv(f, stringsAsFactors = FALSE)
    if (nrow(T) < n_trials) {
      stop("Expected ", n_trials, " rows in ", f, " got ", nrow(T))
    }
    reward_a[i, ] <- T$money_left[seq_len(n_trials)]
    time_a[i, ] <- T$weeks_left[seq_len(n_trials)]
    reward_b[i, ] <- T$money_right[seq_len(n_trials)]
    time_b[i, ] <- T$weeks_right[seq_len(n_trials)]
    ch <- T$choice[seq_len(n_trials)]
    choice[i, ] <- ifelse(ch == "left", 1L, 2L)
  }

  ll <- ifelse(reward_a > reward_b, choice == 1L, choice == 2L)
  participant_order <- order(rowMeans(ll))
  reward_a <- reward_a[participant_order, , drop = FALSE]
  reward_b <- reward_b[participant_order, , drop = FALSE]
  time_a <- time_a[participant_order, , drop = FALSE]
  time_b <- time_b[participant_order, , drop = FALSE]
  choice <- choice[participant_order, , drop = FALSE]

  rLL <- matrix(NA_real_, n_participants, n_trials)
  rSS <- matrix(NA_real_, n_participants, n_trials)
  tLL <- matrix(NA_real_, n_participants, n_trials)
  tSS <- matrix(NA_real_, n_participants, n_trials)
  decision <- matrix(NA_integer_, n_participants, n_trials)
  for (i in seq_len(n_participants)) {
    for (j in seq_len(n_trials)) {
      if (reward_a[i, j] > reward_b[i, j]) {
        decision[i, j] <- as.integer(choice[i, j] == 1L)
        rLL[i, j] <- reward_a[i, j]
        rSS[i, j] <- reward_b[i, j]
        tLL[i, j] <- time_a[i, j]
        tSS[i, j] <- time_b[i, j]
      } else {
        decision[i, j] <- as.integer(choice[i, j] == 2L)
        rLL[i, j] <- reward_b[i, j]
        rSS[i, j] <- reward_a[i, j]
        tLL[i, j] <- time_b[i, j]
        tSS[i, j] <- time_a[i, j]
      }
    }
  }

  list(
    nParticipants = n_participants,
    nTrials = n_trials,
    decision = decision,
    rLL = rLL,
    rSS = rSS,
    tLL = tLL,
    tSS = tSS,
    muPrecScale = mu_prec_scale
  )
}
