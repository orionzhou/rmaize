#' downsample the more abundant class to be balanced with the majority class
#'
#' @export
downsample <- function(ti, seed=1, colname='status') {
  #{{{
  if(colname != 'status')
    ti = ti %>% rename(status = get('colname'))
  levs = levels(ti$status)
  stopifnot(length(levs) == 2)
  tis = ti %>% count(status) %>% arrange(n)
  lev1 = tis %>% pluck('status', 1)
  lev2 = tis %>% pluck('status', 2)
  n1 = tis %>% pluck('n', 1)
  n2 = tis %>% pluck('n', 2)
  ti1 = ti %>% filter(status == lev1)
  set.seed(seed)
  ti2 = ti %>% filter(status == lev2) %>% slice(sample(n2, n1))
  ti1 %>% bind_rows(ti2)
  #}}}
}

#' fit a splited train set and validate in test set
#'
#' @export
fit_split <- function(object, model, split, ...) {
  #{{{
  if (inherits(object, "formula")) {
    object <- add_model(add_formula(workflow(), object, blueprint = hardhat::default_formula_blueprint(indicators = FALSE)), model)
  }
  tune::last_fit(object, split, ...)
  #}}}
}

#' Get tuner
#'
#' @export
get_tuner <- function(alg='rf', cpu=1, mode='classification') {
    #{{{
    if(alg == 'svm') {
      svm_poly(cost=tune(), degree=tune(),
               scale_factor=tune(), margin=tune()) %>%
      set_engine("kernlab") %>%
      set_mode(mode)
    } else if (alg == 'xgb') {
      boost_tree(mtry=tune(), trees=tune(), min_n=tune(),
                 tree_depth=tune(), learn_rate=tune(), loss_reduction=tune()) %>%
      set_engine("xgboost", importance="permutation", nthread=cpu) %>%
      set_mode(mode)
    } else if (alg == 'rf') {
      rand_forest(trees=tune(), min_n=tune()) %>%
      set_engine("ranger", importance="permutation", num.threads=cpu) %>%
      set_mode(mode)
    }
    #}}}
}

#' Get grid
#'
#' @export
get_grid <- function(alg='rf', nlevel=3) {
    #{{{
    if(alg == 'svm') {
    } else if (alg == 'xgb') {
    } else if (alg == 'rf') {
       grid_regular(
                    #finalize(mtry(trans=sqrt_trans()), data),
                    trees(),
                    min_n(),
                    levels = nlevel)
    }
    #}}}
}

#' Get default model
#'
#' @export
get_model <- function(alg='rf', cpu=1, mode='classification') {
    #{{{
    if(alg == 'svm') {
      svm_poly() %>%
      set_engine("kernlab") %>%
      set_mode(mode)
    } else if (alg == 'xgb') {
      boost_tree() %>%
      set_engine("xgboost", importance="permutation", nthread=cpu) %>%
      set_mode(mode)
    } else if (alg == 'rf') {
      rand_forest() %>%
      set_engine("ranger", importance="permutation", num.threads=cpu) %>%
      set_mode(mode)
    }
    #}}}
}

#' Default metrics set for motif prediction
#'
#' @export
metrics6 <- metric_set(sens,spec,precision,accuracy, f_meas, roc_auc, pr_auc)

wfl_fit <- function(wfl) wfl %>% pluck(".workflow", 1) %>% pull_workflow_fit()
wfl_metric <- function(wfl) wfl %>% collect_metrics() %>% select(metric=.metric,estimate=.estimate)
wfl_pred <- function(wfl) wfl %>% pluck('.predictions', 1) %>% rename(pred=.pred_class,truth=status)

#' obtian balanced metrics from truth/prediction
#'
#' @export
metric_balanced <- function(pred) {
  #{{{
  cm = pred %>% conf_mat(truth, pred)
  cm$table = t(t(cm$table) / colSums(cm$table) * 100)
  cm %>% summary() %>% select(metric=.metric, estimate=.estimate)
  #}}}
}

recipe_ml <- function(train_data) {
  #{{{
  rec = recipe(status ~ ., data = train_data) %>%
    #update_role(sid, new_role = "ID") %>%
    #step_medianimpute(all_numeric(), -all_outcomes()) %>%
    #step_modeimpute(all_nominal(), -all_outcomes()) %>%
    #step_novel(all_nominal()) %>%
    step_zv(all_predictors())
    #step_dummy(all_nominal()) %>%
    #themis::step_downsample(all_outcomes())
  rec
  #}}}
}

#' train a simple model using default parameters
#'
#' @export
ml0 <- function(ti, alg='rf', split_prop=.9, cpu=1, seed=26) {
  #{{{
  require(vip)
  set.seed(seed)
  data_split <- initial_split(ti, prop = split_prop)
  train_data <- training(data_split)
  test_data  <- testing(data_split)
  #
  rec = recipe_ml(train_data)
  #
  wfl = workflow() %>%
    add_recipe(rec) %>%
    add_model(get_model(alg=alg, cpu=cpu)) %>%
    fit_split(split = data_split, metrics = metrics6)
  # metrics
  ft = wfl_fit(wfl)
  metric = wfl_metric(wfl)
  pred = wfl_pred(wfl)
  vis = vi(ft) %>% as_tibble()
  #metricB = metric_balanced(pred)
  #
  #tibble(fit=list(ft), vis=list(vis), pred=list(pred), metric=list(metric), metricB=list(metricB))
  tibble(vis=list(vis), pred=list(pred), metric=list(metric))
  #}}}
}

#' train best model using gird search and vfold-CV
#'
#' @export
ml1 <- function(ti, alg='rf', split_prop=.9, fold=10, nlevel=3, cpu=1, seed=26) {
  #{{{
  require(vip)
  set.seed(seed)
  cat("seed =", seed, "\n")
  data_split <- initial_split(ti, prop = split_prop)
  train_data <- training(data_split)
  test_data  <- testing(data_split)
  #
  rec = recipe_ml(train_data)
  #
  wfl = workflow() %>%
    add_recipe(rec) %>%
    add_model(get_tuner(alg=alg, cpu=cpu))
  folds = vfold_cv(train_data, v = fold)
  mfit = wfl %>%
      tune_grid(resamples = folds,
                grid = get_grid(alg=alg, nlevel=nlevel),
                metrics = metric_set(f_meas))
  mres = mfit %>% collect_metrics(summarize = F)
  mfit %>% show_best(metric='f_meas') %>% print(width=Inf)
  param = mfit %>% select_best(metric='f_meas')
  wfl = wfl %>% finalize_workflow(param) %>%
      fit_split(split = data_split, metrics = metrics6)
  # metrics
  ft = wfl_fit(wfl)
  metric = wfl_metric(wfl)
  pred = wfl_pred(wfl)
  vis = vi(ft) %>% as_tibble()
  #metricB = metric_balanced(pred)
  #
  #tibble(fit=list(ft), vis=list(vis), pred=list(pred), metric=list(metric),
         #mres=list(mres), metricB=list(metricB))
  tibble(param=list(param), vis=list(vis), pred=list(pred), metric=list(metric), mres=list(mres))
  #}}}
}
