This repository contains the code for the models. All the models are build in MATLAB (R2015a) and run on a MacBook Pro with a 2.6 GHz quad-core Intel Core i7 processor. In addition several of the models use MATLAB’s Parallel Computing Toolbox with 4 local workers taking advantage of the quad-core processor structure to execute several runs in parallel.

To execute a **single repetition** of the model with specific parameter values (change the values of `pref`) use respectively :

- [ABM_static.m](ABM_static.m)
- [ABM_ind_static.m](ABM_ind_static.m)

To execute the **full experimental design** of the model with a symmetric and unimodal distribution (i.e. **grid sweep method**) use respectively:

- [GS_sticker.m](GS_sticker.m)
- [GS_aggregator.m](GS_aggregator.m)
- [GS_hunter.m](GS_hunter.m)
- [GS_maxcov.m](GS_maxcov.m)
- [GS_maxcovrnd.m](GS_maxcovrnd.m)

To execute the **full experimental design** of the model with an asymmetric and bimodal distribution (i.e. **Monte Carlo parameterisation method**) use respectively:

- [MCP_aggregator.m](MCP_aggregator.m)
- [MCP_hunter.m](MCP_hunter.m)
- [MCP_maxcov.m](MCP_maxcov.m)
- [MCP_maxcovrnd.m](MCP_maxcovrnd.m)
- [MCP_maxcov-inductor.m](MCP_maxcov-inductor.m)
- [MCP_maxcov-inductor-GA.m](MCP_maxcov-inductor-GA.m)