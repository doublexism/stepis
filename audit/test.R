library(rbenchmark)
benchmark(
  test <- test %>% mutate(a = exp(SNP1)),
  test$a <- map_dbl(test$SNP1, exp),
  test$a <- lapply(test$SNP1, exp) %>% c(),
  test$a <- exp(test$SNP1)
)

source("simulation\\simulate data LE.R")
test <-sib_sim(100)

coxph(Surv(age-20, Y) ~ SNP1 + SNP2 + sex  + cov1 + cov2 + cov8, data = test)


mean <- rnorm(1000/2, 0, 0.2) %>% rep(each = 2)
eps <- mvrnorm(500, rep(1,2), diag(rep(0.8, 2))) %>% t() %>% as.vector()
# mean <- mean + rnorm(1000, 0, 0.5)
mean <- mean+eps+1
fid <- rep(1:500,2)
lmer(mean ~(1|fid)) %>% summary()


x <- sample(c(0,1,2), 1000, prob = c(0.25, 0.5, 0.25), replace = TRUE)
y <- sample(c(0,1,2), 1000, prob = c(0.25, 0.5, 0.25), replace = TRUE)
i <- x*y

a <- function(){
  a <- foreach(i=1:3) %do%
     data.frame(sqrt(i)) %>%
  write_rds(a, "data/test.rds")
}
