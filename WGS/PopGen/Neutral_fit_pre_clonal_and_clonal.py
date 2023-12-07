import sys

tumor_id=sys.argv[1]

import math

from pyabc.external.r import R

rname = "~/pyABC/Medulloblastoma/Subclonal_evolution/" + tumor_id + "/Run_model.R"
r = R(rname)

model = r.model("myModel")
distance = r.distance("mySummaryDistance")
observation = r.observation("mySumStatData")


from pyabc import Distribution, RV, ABCSMC, sampler


prior = Distribution(n_clonal=RV("uniform", 0, 10000), mu=RV("uniform", 0.1, 19.9), delta=RV("uniform", 0, 0.99), purity_adj=RV("uniform", -0.05, 0.1))

from pyabc.sampler import MulticoreEvalParallelSampler

sample_specs=sampler.MulticoreParticleParallelSampler(n_procs=8)

abc = ABCSMC(model, prior, distance, sampler = sample_specs, population_size = 1000)

import os
from tempfile import gettempdir

db = "sqlite:///" + "/home/koerber/pyABC/Medulloblastoma/Subclonal_evolution/" + tumor_id + "/" + tumor_id + ".db"

abc.new(db, r.observation("mySumStatData"))

history = abc.run(minimum_epsilon=1, max_nr_populations=25)


