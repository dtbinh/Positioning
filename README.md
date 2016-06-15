## Competitive location behaviour with foresight
Jonas K. Sekamane *(Supervised by Johan Lagerlöf)*.

**Abstract:** This paper considers a competitive location model, where firms compete with one another through their location decision in a two-dimensional market space. The paper studies how firms locate and the  corresponding competitive environment. In particular the paper investigates when firms agglomerate at the centre or seek out niche segments of the market. And whether one firm captures a predominant share of customers, when firms compete through location, or if the market is evenly partitioned among firms. Both a symmetric unimodal and an asymmetric and bimodal distribution of consumers are examined. An agent-based modelling approach is used and the paper proposes a new decision rule with foresight -- as an alternative to the current heuristic decision rules, where firms are oblivious to the simultaneous action and future response from competing firms. The new decision rule uses *inductive reasoning*, i.e. the firm holds several hypotheses on how competing firms chooses to locate, and use the most probable hypothesis to predict the future location of competing firms. In addition the firm gradually discards poorly performing hypotheses and forms new hypotheses, thus the hypotheses of the firms are endogenous to the model. Firms acting with foresight enable strategic considerations.

[Full paper](https://github.com/jsekamane/Positioning/raw/master/Paper/sekamane_location.pdf)

* * *

This repository contains the code for the models. All the models are build in MATLAB (R2015a) and run on a MacBook Pro with a 2.6 GHz quad-core Intel Core i7 processor. In addition several of the models use MATLAB’s Parallel Computing Toolbox with 4 local workers taking advantage of the quad-core processor structure to execute several runs in parallel.

To execute a **single repetition** of the model with specific parameter values (change the values of `pref`) use respectively :

- [ABM_static.m](/Models/ABM/ABM_static.m)
- [ABM_ind_static.m](/Models/ABM/ABM_ind_static.m)

To execute the **full experimental design** of the model with a symmetric and unimodal distribution (i.e. **grid sweep method**) use respectively:

- [GS_sticker.m](/Models/ABM/GS_sticker.m)
- [GS_aggregator.m](/Models/ABM/GS_aggregator.m)
- [GS_hunter.m](/Models/ABM/GS_hunter.m)
- [GS_maxcov.m](/Models/ABM/GS_maxcov.m)
- [GS_maxcovrnd.m](/Models/ABM/GS_maxcovrnd.m)

To execute the **full experimental design** of the model with an asymmetric and bimodal distribution (i.e. **Monte Carlo parameterisation method**) use respectively:

- [MCP_aggregator.m](/Models/ABM/MCP_aggregator.m)
- [MCP_hunter.m](/Models/ABM/MCP_hunter.m)
- [MCP_maxcov.m](/Models/ABM/MCP_maxcov.m)
- [MCP_maxcovrnd.m](/Models/ABM/MCP_maxcovrnd.m)
- [MCP_maxcov-inductor.m](/Models/ABM/MCP_maxcov-inductor.m)
- [MCP_maxcov-inductor-GA.m](/Models/ABM/MCP_maxcov-inductor-GA.m)