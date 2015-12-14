summary.climate_histories = function(object,dim=1,probs=c(0.05,0.95),...) {

# Output the climate histories for time grid point

out = apply(object$histories[,,dim],2,'quantile',prob=probs)
df = data.frame('time_grid'=object$time_grid,t(out))
print(df)

invisible(df)

}
