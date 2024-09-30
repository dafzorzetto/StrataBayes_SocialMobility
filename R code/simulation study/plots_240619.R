OUT_1=model_setting_2_correct[[1]]
scenario=scenario_2

par(mfrow=c(3,1))
for(i in 1:3){
  plot(OUT_1$BETA_0_POST[i,], type="l", ylab = paste0("Beta_0_Post ", i))
  abline(h=mean(OUT_1$BETA_0_POST[i,]), col="red")
  abline(h=scenario[[1]]$parameters$beta_0[i], col="blue")
}

par(mfrow=c(3,1))
for(i in 1:3){
  plot(OUT_1$BETA_1_POST[i,], type="l", ylab = paste0("Beta_1_Post ", i))
  abline(h=mean(OUT_1$BETA_1_POST[i,]), col="red")
  abline(h=scenario[[1]]$parameters$beta_1[i], col="blue")
}

par(mfrow=c(2,1))
plot(OUT_1$SIGMA_0_POST, type="l")
abline(h=mean(OUT_1$SIGMA_0_POST), col="red")
abline(h=scenario[[1]]$parameters$sigma_0, col="blue")
plot(OUT_1$SIGMA_1_POST, type="l")
abline(h=mean(OUT_1$SIGMA_1_POST), col="red")
abline(h=scenario[[1]]$parameters$sigma_1, col="blue")

par(mfrow=c(3,4))
for(v in 1:3){
  for(i in 1:4){
    plot(OUT_1$ETA_0_POST[v,i,], type="l",
         main=paste0("cluster ", v, " eta_0_",i))
    abline(h=mean(OUT_1$ETA_0_POST[v,i,]), col="red")
    abline(h=scenario[[1]]$parameters$eta_0[v,i], col="blue")
  }
}


par(mfrow=c(3,5))
for(v in 1:3){
  for(i in 1:5){
  plot(OUT_1$ETA_1_POST[v,i,], type="l",
       main=paste0("cluster ", v, " eta_1_",i))
  abline(h=mean(OUT_1$ETA_1_POST[v,i,]), col="red")
  abline(h=scenario[[1]]$parameters$eta_1[v,i], col="blue")
  }
}

par(mfrow=c(3,1))
for(i in 1:3){
  plot(OUT_1$SIGMA_Y_0_POST[i,], type="l")
  abline(h=mean(OUT_1$SIGMA_Y_0_POST[i,]), col="red")
  abline(h=scenario[[1]]$parameters$sigma_y_0[i], col="blue")
}

par(mfrow=c(3,1))
for(i in 1:3){
  plot(OUT_1$SIGMA_Y_1_POST[i,], type="l")
  abline(h=mean(OUT_1$SIGMA_Y_1_POST[i,]), col="red")
  abline(h=scenario[[1]]$parameters$sigma_y_1[i], col="blue")
}


#scenario 2
par(mfrow=c(2,2))
hist(scenario[[1]]$simulated_full$P_0, breaks=100, xlim=c(-5,3), main="P0 Simulated", xlab=NULL, ylab=NULL)
hist(scenario[[1]]$simulated_full$P_1, breaks=100, xlim=c(-7,6), main="P1 Simulated", xlab=NULL, ylab=NULL)
hist(rowMeans(OUT_1$P0_POST), breaks=100, xlim=c(-5,3), main="P0 Post", xlab=NULL, ylab=NULL)
hist(rowMeans(OUT_1$P1_POST), breaks=100, xlim=c(-7,6), main="P1 Post", xlab=NULL, ylab=NULL)

par(mfrow=c(2,2))
hist(scenario[[1]]$simulated_full$Y_0, breaks=100, xlim=c(-15,15), main="Y0 Simulated", xlab=NULL, ylab=NULL)
hist(scenario[[1]]$simulated_full$Y_1, breaks=100, xlim=c(-10,15), main="Y1 Simulated", xlab=NULL, ylab=NULL)
hist(rowMeans(OUT_1$Y0_POST), breaks=100, xlim=c(-15,15), main="Y0 Post", xlab=NULL, ylab=NULL)
hist(rowMeans(OUT_1$Y1_POST), breaks=100, xlim=c(-10,15), main="Y1 Post", xlab=NULL, ylab=NULL)
