import Banmi
import Mix

x = readcsv("stlouis.csv")

desc = Mix.describe(x, 3)
theta = Mix.em(desc)
genloc = Mix.getparam(desc, theta)

Mix.rngseed(1234)
wimp, zimp = Mix.impute(desc, theta, x)

# desc = Banmi.describe(x, 3)
# dp_weight = 10.0
# lambda_a = 3.0
# lambda_b = 2.0
# hp = Banmi.hyperparameters(desc, genloc, dp_weight, lambda_a, lambda_b)

println(wimp, zimp)
