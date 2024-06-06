**HEAT ECHANGER NETWORKS**

The heat exchanger network (HEN) problem focuses on improving energy savings in industrial sites. HENs aim to integrate the heat from hot and cold process streams to reduce reliance on external energy sources, such as heating and cooling utilities, thereby increasing the energy efficiency of the process. A proper balance must be achieved between energy savings, the number of heat exchangers, the required area to recover that energy, and the network structure to arrive at the final design with the best possible total annual cost (TAC). However, since HENs naturally result in a non-convex, non-linear, combinatorial optimization problem, developing low-cost feasible HENs is particularly challenging. For these reasons, HENs are one of the most important problems in process system engineering and have been an open subject of investigation for the last five decades.

Originally, this problem was approached using a thermodynamic perspective and elements of graph theory to develop the well-known pinch-point methodology. Despite being widely used in the industry for their simplicity, these approaches require a large number of manual calculations. This approach was systematized through the use of sequential mathematical programming models solved individually, originally proposed as a linear programming (LP) model to determine only the maximum amount of recoverable energy, then a mixed-integer linear programming (MILP) model to establish the minimum number of matches between process streams, and a nonlinear programming (NLP) model to determine the TAC.

**Approach Based on Mathematical Programming Using Transportation Models**

One way to formulate the problem involves decomposing the temperature range, considering only the input and output temperatures of all streams, to establish a temperature scale whose boundaries define intervals. The hot temperature intervals serve as supply nodes in the transportation model, where all the available energy in these nodes must be exchanged or transported to the demand nodes represented by the cold temperature intervals. According to the first law of thermodynamics, each hot process stream can exchange heat from a given temperature interval to a stream in another cold temperature interval that is at an equal or lower temperature level. The position of the intervals depends on the input and output temperatures of the process streams and the heat recovery approach temperature (HRAT) parameter, which defines the relative displacement between the temperature scales of the hot and cold process streams. A low HRAT parameter value implies a greater amount of vertical heat exchange, translating to an increase in heat recovery capacity, but at the cost of smaller temperature differences, leading to a higher heat exchange area requirement. Conversely, as the HRAT parameter value increases, the opposite occurs. Fig. 1 shows a simplified diagram representing two contiguous temperature intervals for two hot and two cold process streams.

The proposed novel models address this problem simultaneously, tackling all elements involved in the heat recovery process in one step. These approaches aim to achieve a proper balance between the amount of recovered heat, the number of matches, the heat exchange area, and the TAC.

Fig. 1. Feasible thermal energy flows in two contiguous temperature intervals.

**Model (M) – MILP**

(1)

Subject to

*Supply and demand constraints*

(2)

(3)

*Logical Constraints*

(4)

(5)

(6)

Where M*i,j*, Mhu*j*, Mcu*i* is the maximum amount of heat that can be exchanged between two process streams or utilities.

*Non-negativity constraints*

(7)

(8)

**Model (N) – NLP**

(1)

where Nhx is the user-defined number of heat recovery units..

Subject to

*Supply and demand constraints*

(2)

(3)

*Heat exchange area*

(4)

(5)

(6)

*Number of heaters and coolers*

(7)

(8)

Where ε is a very small number

*Non-negativity constraints*

(9)

(10)

**Model (L) - LP**

(1)

*Supply and demand constraints*

(2)

(3)

Cotas para la transferencia de calor

(4)

(5)

(6)

*Special linearized coefficients for model (L)*

**Results**

All calculations were performed in a personal DELL XPS 8500 CPU Intel i7-3770 3.7 GHz with 16 GB of RAM memory running windows. The models were implemented and solved using general algebraic modeling system (GAMS), using CPLEX as LP solver, CONOPT as NLP solver, and DICOPT as MILP solver.

Based on the process stream and cost information, the previously shown models allow for the determination of the TAC for a given HRAT. This can be systematized to determine the optimal value of the HRAT parameter that provides the lowest TAC. Table 1 shows an illustrative example of a problem with 11 hot process streams and two cold streams, as well as information related to auxiliary services and costs.

*Table 1. Process streams and costs information for the case study addressed.*

| Stream                                            | TS (K) | TO (K) | F (kW K-1) | h (kW m-2 K-1) |
|---------------------------------------------------|--------|--------|------------|----------------|
| H1                                                | 523.75 | 363.15 | 132.20     | 0.26           |
| H2                                                | 413.35 | 312.65 | 106.50     | 0.26           |
| H3                                                | 576.75 | 543.35 | 234.98     | 0.41           |
| H4                                                | 563.15 | 388.15 | 39.81      | 0.47           |
| H5                                                | 483.15 | 436.15 | 115.76     | 0.33           |
| H6                                                | 521.95 | 383.15 | 31.81      | 0.72           |
| H7                                                | 550.15 | 395.05 | 24.58      | 0.57           |
| H8                                                | 443.25 | 333.15 | 33.93      | 0.45           |
| H9                                                | 451.75 | 382.05 | 47.85      | 0.6            |
| H10                                               | 633.15 | 563.15 | 39.81      | 0.47           |
| H11                                               | 632.75 | 553.15 | 24.53      | 0.47           |
| C1                                                | 403.15 | 623.15 | 289.92     | 0.72           |
| C2                                                | 303.15 | 403.15 | 202.48     | 0.26           |
| HU                                                | 773.15 | 772.15 | –          | 0.53           |
| CU                                                | 293.15 | 313.15 | –          | 0.53           |
| Heat exchanger cost (\$ yr-1) = 2500+55[Area(m2)] |        |        |            |                |
| Hot utility cost = 100 (\$ kW-1 yr-1)             |        |        |            |                |
| Cold utility cost = 10 (\$ kW-1 yr-1)             |        |        |            |                |

Fig 2. TAC for different selected HRAT values.

**Conclusions**

The proposed models allow for the estimation of the TAC of a heat exchanger network problem in a single step, considering only the process stream data and costs, for any given HRAT value. Unlike traditional methods such as supertargeting, the proposed models do not require manual calculations. Particularly, the LP (Linear Programming) model is a promising alternative for this type of problem because formulations based on mixed-integer programming are characterized as strongly NP-hard problems. This implies that it is unlikely to find an algorithm that solves these problems in polynomial computational time, leading to a combinatorial explosion in larger-scale problems.

**Nomenclature**

*Indices*

|   | Index for a cold process stream                              |
|---|--------------------------------------------------------------|
|   | Index for a hot process stream                               |
|   | Index for hot temperature intervals and interval boundaries  |
|   | Index for cold temperature intervals and interval boundaries |

*Sets*

|      | Hot streams                            |
|------|----------------------------------------|
|      | Cold Streams                           |
|      | Hot temperature intervals {1,2,…,NOL}  |
| *L1* | Cold temperature intervals {1,2,…,NOL} |

*Parameters*

|   | Heat recovery approach temperature (°C)                                                              |
|---|------------------------------------------------------------------------------------------------------|
|   | Supply temperature (°C)                                                                              |
|   | Target temperature (°C)                                                                              |
|   | Supply temperature for hot and cold utility (°C)                                                     |
|   | Target temperature for hot and cold utility (°C)                                                     |
|   | Hot process stream temperature on the upper boundary of interval (°C)                                |
|   | Cold process stream temperature on the upper boundary of interval (°C)                               |
|   | Heat capacity flow rate (kW °C-1)                                                                    |
|   | Heat film transfer coefficient (kW °C-1 m-2)                                                         |
|   | Heat available from hot stream  in interval  (kW)                                                    |
|   | Heat demand of cold stream in interval *m* (kW)                                                      |
|   | Total number of temperature intervals                                                                |
|   | Penalty parameter (e.g. )                                                                            |
|   | Transportation cost of heat recovery                                                                 |
|   | Transportation cost of heating utilities                                                             |
|   | Transportation cost of cooling utilities                                                             |
|   | Upper bound (Big-M parameter) for heat exchange between hot stream *I* and cold stream *j* (kW)      |
|   | Upper bound (Big-M parameter) for hot utility requirements of stream *j* (kW)                        |
|   | Upper bound (Big-M parameter) for cold utility requirements of stream *i* (kW)                       |
|   | Logarithmic mean temperature difference between temperature intervals involving process streams (°C) |
|   | Logarithmic mean temperature difference for cold temperature intervals involving hot utilities (°C)  |
|   | Logarithmic mean temperature difference for hot temperature intervals involving cold utilities (°C)  |
|   | Fixed charge of capital cost (\$)                                                                    |
|   | Capital cost law coefficient (\$m−2βk)                                                               |
|   | Capital cost law exponent                                                                            |
|   | Number of heat exchangers                                                                            |

*Positive continuous variables*

|   | Heat exchanged between hot stream  in the interval *l* to the cold stream in the interval *m* (kW) |
|---|----------------------------------------------------------------------------------------------------|
|   | Heat utility required of a cold streamin the interval *m* (kW)                                     |
|   | Cold utility required of a hot stream *i* in the interval *l* (kW)                                 |

*Binary variables*

|   | Binary variable for heat match between hot stream  and cold stream    |
|---|-----------------------------------------------------------------------|
|   | Binary variable for heat match between hot utilities and cold stream  |
|   | Binary variable for heat match between hot stream  and cold utilities |

**Calculable parameters**

*Enthalpies of hot process streams*

*Enthalpies of cold process streams*

*Logarithmic mean temperature difference*

\\
