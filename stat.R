library(stringr)
library(ROCR)
library(ggplot2)
library(PropCIs)

### 01 find optimal cutoff ####
find_optimal_cutoff = function(prob, label) {
  val_pred = ROCR::prediction(prob, label)
  val_perf = ROCR::performance(val_pred, 'tpr', 'fpr')
  val_table = data.frame(fpr=unlist(val_perf@x.values), tpr=unlist(val_perf@y.values),
                         cutoff=unlist(val_perf@alpha.values))
  val_table[,'youden'] = (1-val_table[,'fpr'])+val_table[,'tpr']
  optimal_cutoff = val_table[which.max(val_table[,'youden']), 'cutoff']
  
  pred_cutoff = ifelse(prob >= optimal_cutoff, 1, 0)
  cf_matrix = table(pred_cutoff, label)
  sens = cf_matrix[2,2] / sum(cf_matrix[,2])
  spec = cf_matrix[1,1] / sum(cf_matrix[,1])
  ppv = cf_matrix[2,2] / sum(cf_matrix[2,])
  npv = cf_matrix[1,1] / sum(cf_matrix[1,])
  return(list(optimal_cutoff=optimal_cutoff, cf_matrix=cf_matrix,
              sens=sens, spec=spec, ppv=ppv, npv=npv))
}


### 02 find cutoff point based on sensitivity or specificity ####
find_cutoff_point = function(prob, label, type, value) {
  if (type == 'tpr') {other_type = 'tnr'}
  else if (type == 'tnr') {other_type = 'tpr'}
  
  val_pred = ROCR::prediction(prob, label)
  val_perf = ROCR::performance(val_pred, 'tpr', 'fpr')
  val_table = data.frame(tnr=1-unlist(val_perf@x.values), tpr=unlist(val_perf@y.values),
                         cutoff=unlist(val_perf@alpha.values))
  val_table[,'youden'] = val_table[,'tnr']+val_table[,'tpr']
  
  val_table[,'point'] = abs(val_table[,type] - value)
  points_df = val_table[val_table[,'point'] == min(val_table[,'point']),]
  cutoff_point = points_df[which.max(points_df[,other_type]), 'cutoff']
  
  pred_cutoff = ifelse(prob >= cutoff_point, 1, 0)
  cf_matrix = table(pred_cutoff, label)
  sens = cf_matrix[2,2] / sum(cf_matrix[,2])
  spec = cf_matrix[1,1] / sum(cf_matrix[,1])
  ppv = cf_matrix[2,2] / sum(cf_matrix[2,])
  npv = cf_matrix[1,1] / sum(cf_matrix[1,])
  return(list(cutoff_point=cutoff_point, type=type, value=value, cf_matrix=cf_matrix,
              sens=sens, spec=spec, ppv=ppv, npv=npv))
}


### 03 Calculate diagnostic performance ####
find_cf_stat = function(prob, label, cutoff) {
  pred_cutoff = ifelse(prob >= cutoff, 1, 0)
  cf_matrix = table(pred_cutoff, label)
  
  pos_d = cf_matrix[2,2]
  neg_d = cf_matrix[1,1]
  
  pos_n = cf_matrix[,2]
  neg_n = cf_matrix[,1]
  
  pos_pred = cf_matrix[2,]
  neg_pred = cf_matrix[1,] 
  
  sens_value = pos_d/sum(pos_n)
  spec_value = neg_d/sum(neg_n)
  
  sens_value = cf_matrix[2,2]/sum(cf_matrix[,2])
  spec_value = cf_matrix[1,1]/sum(cf_matrix[,1])
  
  acc_nume = sum(diag(cf_matrix))
  acc_denomi = sum(cf_matrix)
 
  sens = paste0(sprintf('%.1f', 100*sens_value), '%')
  sens_count = paste0(as.character(pos_d), '/', sum(pos_n))
  sens_ci = stringr::str_c(sprintf('%.1f', 100*PropCIs::exactci(pos_d, sum(pos_n), conf.level=0.95)$conf.int), '%')
  
  spec = paste0(sprintf('%.1f', 100*spec_value), '%')
  spec_count = paste0(as.character(neg_d), '/', sum(neg_n))
  spec_ci = stringr::str_c(sprintf('%.1f', 100*PropCIs::exactci(neg_d, sum(neg_n), conf.level=0.95)$conf.int), '%')
  
  ppv = paste0(sprintf('%.1f', 100*pos_d/sum(pos_pred)), '%')
  ppv_count = paste0(as.character(pos_d), '/', sum(pos_pred))
  ppv_ci = stringr::str_c(sprintf('%.1f', 100*PropCIs::exactci(pos_d, sum(pos_pred), conf.level=0.95)$conf.int), '%')
  
  npv = paste0(sprintf('%.1f', 100*neg_d/sum(neg_pred)), '%')
  npv_count = paste0(as.character(neg_d), '/', sum(neg_pred))
  npv_ci = stringr::str_c(sprintf('%.1f', 100*PropCIs::exactci(neg_d, sum(neg_pred), conf.level=0.95)$conf.int), '%')
  

  plr = sens_value / (1 - spec_value)
  plr_value = sprintf('%.2f', sens_value / (1 - spec_value))
  plr_se = sqrt(1/pos_d - 1/sum(pos_n) + 1/cf_matrix[2,1] - 1/sum(neg_n))
  plr_ci = stringr::str_c(sprintf('%.2f', exp(log(plr)+(c(-1,1)*qnorm(0.975)*plr_se))))
  
  nlr = (1-sens_value) / (spec_value)
  nlr_value = sprintf('%.2f', (1-sens_value) / (spec_value))
  nlr_se = sqrt(1/cf_matrix[1,2] - 1/sum(pos_n) + 1/neg_d - 1/sum(neg_n))
  nlr_ci = stringr::str_c(sprintf('%.2f', exp(log(nlr)+(c(-1,1)*qnorm(0.975)*nlr_se))))
  
  output_df = data.frame(value=c(sens, spec, ppv, npv, plr_value, nlr_value),
                         count=c(sens_count, spec_count, ppv_count, npv_count, NA, NA),
                         ci_lower = c(sens_ci[1], spec_ci[1], ppv_ci[1], npv_ci[1], plr_ci[1], nlr_ci[1]),
                         ci_upper = c(sens_ci[2], spec_ci[2], ppv_ci[2], npv_ci[2], plr_ci[2], nlr_ci[2]))
  rownames(output_df) = c('Sensitivity', 'Specificity', 'PPV', 'NPV', 'PLR', 'NLR')
  return(output_df) 
}


### 04 Draw multiple ROC curve ####
compare_roc_curve_ci = function(xlsx_list, prob_name, label_name, legend_color, legend_lty,
                                model_info, title_name, legend_key_width, png_name) {
                            
  pred_list = list()
  perf_list = list()
  auc_list = list()
 
  for (i in 1:length(xlsx_list)) {
    pred_list[[i]] = ROCR::prediction(xlsx_list[[i]][,prob_name], xlsx_list[[i]][,label_name])
    perf_list[[i]] = ROCR::performance(pred_list[[i]], 'tpr', 'fpr')
    auc_list[[i]] = sprintf('%.3f', pROC::ci.auc(xlsx_list[[i]][,label_name], xlsx_list[[i]][,prob_name], method='delong'))
  }
  
  roc_curve_theme = ggplot2::theme(
    title = ggplot2::element_text(size=12),
    axis.title.x = ggplot2::element_text(size=10),
    axis.title.y = ggplot2::element_text(size=10),
    axis.text.x = ggplot2::element_text(size=8),
    axis.text.y = ggplot2::element_text(size=8),
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.margin = unit(c(0.2,0.2,0.2,0.2), 'cm'),
    legend.position = c(1,0),
    legend.justification = c(1,0),
    legend.text = ggplot2::element_text(size=9.6),
    legend.title = ggplot2::element_text(size=9.6),
    legend.background = ggplot2::element_rect(fill="gray80", color="black"),
    legend.margin = margin(0.1, 0.1, 0.1, 0.1, 'cm'),
    legend.spacing = ggplot2::unit(0.001, "line"),
    legend.key.height = ggplot2::unit(1.5, "line"),
    legend.key.width = ggplot2::unit(1.5, "line")
    # legend.key.width = ggplot2::unit(legend_key_width, 'line')
  )
  
  r_curve = ggplot2::ggplot(data=NULL) + 
    ggplot2::coord_equal(expand=FALSE, xlim=c(-0.01, 1), ylim=c(0, 1.01))+
    ggplot2::xlab('1-Specificity')+
    ggplot2::ylab('Sensitivity')+
    ggplot2::ggtitle(title_name)+
    ggplot2::theme_bw() +
    roc_curve_theme
  
  legend_text = c()
  
  for (j in 1:length(xlsx_list)) {
    each_curve = ggplot2::geom_line(ggplot2::aes_string(x=unlist(perf_list[[j]]@x.values), y=unlist(perf_list[[j]]@y.values), 
                                                        color=as.factor(j), lty=as.factor(j)),  size=1)
    r_curve = r_curve + each_curve
    legend_text[j] = stringr::str_c(model_info[j], ': ', auc_list[[j]][2], ' ', '(', auc_list[[j]][1], '-', auc_list[[j]][3], ')')
  }
  r_curve = r_curve + ggplot2::scale_color_manual(name='AUC (95% C.I.)', values=legend_color, labels=legend_text)
  r_curve = r_curve + ggplot2::scale_linetype_manual(name='AUC (95% C.I.)', values=legend_lty, labels=legend_text)
  r_curve = r_curve + ggplot2::scale_x_continuous(breaks=seq(0, 1, 0.2))
  r_curve = r_curve + ggplot2::scale_y_continuous(breaks=seq(0, 1, 0.2))
  
  print(legend_text)
  return(r_curve)
}
