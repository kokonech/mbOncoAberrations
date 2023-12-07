import math

from pyabc.external import R

r = R("~/pyABC/Medulloblastoma/Expansion_decay_continuous_evol.R")

model = r.model("myModel")
distance = r.distance("mySummaryDistance")
sum_stat = r.summary_statistics("mySummaryStatistics")

observation = r.observation("mySumStatData")


from pyabc import Distribution, RV, ABCSMC

prior = Distribution(delta1=RV("uniform", 0, 0.9), N=RV("uniform", 3, 6), delta2=RV("uniform", 1.0001, 0.4999), psurv=RV("uniform", 0.01, 0.49), muD1=RV("uniform", -9, 5), muD2 = RV("uniform", -9, 5), mu = RV("uniform", 1, 14), r = RV("uniform", 0, 0.999))

from pyabc.sampler import MulticoreEvalParallelSampler

multi_sampler = MulticoreEvalParallelSampler(n_procs=16)

abc = ABCSMC(model, prior, distance, summary_statistics = sum_stat, sampler = multi_sampler, population_size = 1000)

import os
from tempfile import gettempdir

db = "~/pyABC/Medulloblastoma/Expansion_decay_continuous_evol.db"

abc.new(db, r.observation("mySumStatData"))

history = abc.run(minimum_epsilon=0.05, max_nr_populations=25, n_procs=16)
