#status/incomplete #topic/biosensing #topic/SPR #type/index 

date of creation: 2026-05-27

- [ ] Make sure that all assumptions are listed. 

# Models, In Summary
Shown as implemented in the application. Details in the following section
## 1. First-order, one-to-one
$$R(t) = \frac{R_{max} C k_{a}}{Q}\left(1 - e^{-Qt}\right) + s_0e^{-Qt}$$
Where $Q = k_{a}C + k_{d}$  
## 2. Mass-transport Limited
$$\begin{align}
\frac{dR}{dt} &= k_a C_s(R_{max} - R) - k_dR \\
C_s &= \frac{k_{tr}C + k_d R}{k_{tr} + k_a(R_{max} - R)}
\end{align}$$
## 3. Heterogeneous Ligand (one-to-two)
$$\begin{align}
\frac{dR_1}{dt} &= k_{a,1}C(R_{max,1} - R_1) - k_{d,1}R_1\\
\frac{dR_2}{dt} &= k_{a,2}C(R_{max,2} - R_2) - k_{d,2}R_2 \\
R_{obs} &= R_1 + R_2
\end{align}$$
## 4. Conformational Change
$$\begin{align}
\frac{dR_1}{dt} &= k_{a,1}C(R_{max}-R_1 - R_2) - k_{d,1}R_1 - k_{a,2}R_1 + k_{d,2}R_2 \\
\frac{dR_2}{dt} &= k_{a,2}R_1 - k_{d,2}R_2 \\
R &= R_1 + R_2
\end{align}$$
## 5. Bivalent Analyte
$$\begin{align}
\frac{dR_1}{dt} &= 2k_{a,1}C(R_{max} - R_1 - 2R_2) - k_{d,1}R_1 - k_{a,2}R_1(R_{max} - R_1 - 2R_2) + 2k_{d,2}R_2 \\
\frac{dR_2}{dt} &= 2k_{a,2}R_1(R_{max} - R_1 - 2R_2) - 2k_{d,2}R_2 \\
R&=R_1+R_2
\end{align}$$
# Models, In Detail

## Regarding the relationship between response and mass of bound complex:
The 1:1 Langmuir model is linear in surface species, so the response variable $R$ (in RU) is related to the underlying surface concentration $s$ by a constant calibration factor $\alpha$ that absorbs into $R_{max}$​. The ODE retains its form: $$\frac{dR}{dt} = k_{a} C (R_{max} - R) - k_{d} R$$Fitted rate constants $k_{on}$​ and $k_{off}$​ carry their physical units (M$^{−1}$s$^{−1}$ and s$^{-1}$); 
$R_{max}$​ is instrument/surface-specific. This property is preserved for all bread-and-butter models *except* the bivalent analyte model, where the bilinear surface term causes $k_{a2}$​ to inherit calibration-dependent units of RU$^{-1}$ s$^{-1}$.
## Common Variable and Parameter Definitions
The following will be common to all models, unless otherwise specified.

| Table A   | Meaning                     | Unit           |
| --------- | --------------------------- | -------------- |
| $R$       | instrument response         | $RU$           |
| $R_{max}$ | maximum instrument response | $RU$           |
| $s$       | mass of bound complex       |                |
| $C$       | concentration of analyte    | $M$            |
| $k_a$     | association rate constant   | $M^{-1}s^{-1}$ |
| $k_d$     | dissociation rate constant  | $s^{-1}$       |
| $t$       | elapsed time                | seconds        |

Note: $R{max}$ is specific to the instrument, the surface density, and the analyte. It should not be assumed to have physical relevance. It must be a fitted parameter. 
## 1. First-order, one-to-one
#### Variables, Parameters, and Units
All in agreement with Table A.
#### Differential Equation

$$\frac{dR}{dt} = k_{a}C(R_{max} - R) - k_{d}R$$
#### Integrated 
$$R(t) = \frac{R_{max} C k_{a}}{Q}\left(1 - e^{-Qt}\right) + s_0e^{-Qt}$$
Where $Q = k_{a}C + k_{d}$  
## 2. Two Compartment Model of Mass Transport-Limited Binding
The two compartment model is only a good approximation of *mild* to *moderate* mass-transport limited kinetics. Not applicable when the microscopic transport step causes significant concentration gradients in sensing volume, ie when it is not longer reasonable to define the sensing volume as well-mixed.
#### Variables, Parameters, and Units
| Table B  | Meaning                                           | Unit                   |
| -------- | ------------------------------------------------- | ---------------------- |
| $C_{s}$  | concentration of analyte near the binding surface | $M$                    |
| $k_{tr}$ | transport rate constant                           | $RU\cdot M^{-1}s^{-1}$ |
Note that a *lower* value of $k_{tr}$ indicates a $greater$ degree of mass-transport limitation. 
#### Differential Equation
$$\begin{align}
\frac{dC_{s}}{dt} &= k_{tr}(C - C_{s}) - \frac{ds}{dt} \\
\frac{ds}{dt} &= k_{a}C_{s}(s_{max} - s) - k_{d}s
\end{align}$$
Note: This equation remains in $s$ so that the units remain consistent. In this formulation, the unit on $k_{tr}$ is $s^{-1}$. Units of $s$ are also $M$. The transition from $s$ to $R$ gives us the different units of $k_{tr}$, as it absorbs the proportionality.   
#### Quasi-Steady-State Approximation
$C_s$ is expected to reach steady-state far sooner than R, and therefore we can approximate $\frac{dC_s}{dt}=0$, and as a result, we may arrive at the following, simplified form: 
$$\begin{align}
\frac{dR}{dt} &= k_a C_s(R_{max} - R) - k_dR \\
C_s &= \frac{k_{tr}C + k_d R}{k_{tr} + k_a(R_{max} - R)}
\end{align}$$
This approximation is what we use in the data simulator. 
## 3. Heterogeneous Ligand
Models two independent populations of binding sites. This has been our "one-to-two" model. Binding to one type of site is assumed to have no effect whatsoever on the dynamics between the analyte and the other class of binding site. Therefore, the observed signal is the sum of two, one-two-one dynamics. 
#### Variables, Parameters, and Units
In this model $k_{a,1}$ and $k_{a,2}$ share the same unit as $k_a$ in Table A. Likewise, $k_{d,1}$ and $k_{d.2}$ share the same unit as $k_d$ in the same. 
#### Differential Equation
$$\begin{align}
\frac{dR_1}{dt} &= k_{a,1}C(R_{max,1} - R_1) - k_{d,1}R_1\\
\frac{dR_2}{dt} &= k_{a,2}C(R_{max,2} - R_2) - k_{d,2}R_2 \\
R_{obs} &= R_1 + R_2
\end{align}$$
Note that the data simulator employs the integrated form of these equations, as in equation 1.
## 4. Conformational Change
This model represents a single receptor that, upon binding, adopts a second conformation. The second conformation is something like the two species "settling in" to a lower-energy configuration. The second conformation "locks" the two together, and are not modeled to separate until the the pair reverts to the original conformation. Put another way, upon initial binding one of two things can happen: either the analyte releases back into solution, or the two may enter the second conformation. 

Practically, this results in an apparently slow dissociation phase.
#### Variables, Parameters, and Units
| Table C   | Meaning                                         | Unit     |
| --------- | ----------------------------------------------- | -------- |
| $k_{a,2}$ | rate parameter for adopting second conformation | $s^{-1}$ |
Note: The rate at which the complex adopts the second configuration is not dependent on concentration.
#### Differential Equation
- [ ] Standardize $s$ to $R$ #BS2026 
$$\begin{align}
\frac{ds_1}{dt} &= k_{a,1}A(s_{max}-s_1 - s_2) - k_{d,1}s_1 - k_{a,2}s_1 + k_{d,2}s_2 \\
\frac{ds_2}{dt} &= k_{a,2}s_1 - k_{d,2}s_2 \\
s &= s_1 + s_2
\end{align}$$
### Regarding Calculation of "Apparent KD"
$$K_{D,\text{app}}=\frac{k_{d,1}}{k_{a,1}}\cdot \frac{k_{d,2}}{k_{d,2}+k_{a,2}}$$

Factor one: Dissociation constant of initial binding step alone. 
Factor two: fraction of complex that remains in the dissociable AB state at equilibrium
## 5. Bivalent Analyte - 
This model illustrates the case where one species may bind to two receptors at the same time. This requires that the binding species have two distinct appendages or binding regions. This may describe two arms of a y-shaped antibody that can each bind an antigen, or the two binding regions of enzyme that binds a membrane-bound, receptor tyrosine kinase. 

In any case, at any time the total mass of bound receptor, analyte complex may be comprised of some amount of singly-bound or doubly-bound analyte. The probability of a second binding event is greatly dependent on the density of receptor, ie there must second binding site sufficiently close to the first. 

Note that this set of differential equations is not linear in mass of bound complex, therefore there is not a clean linear relationship between mass of bound complex and instrument response. This is absorbed by the second association rate constant, which has units in terms of $RU$. Practically, this means that $k_{a,2}$, much like $R_{max}$ above, depends on experimental conditions, and cannot be assumed to have physical meaning. 

In general, this model has documented issues with parameter identifiability and experimental reproducibility. This will be further explored in a future update. 
- [ ] Remember to expand on this #BS2026
#### Variables, Parameters, and Units
| Table D   | Meaning                                       | Unit            |
| --------- | --------------------------------------------- | --------------- |
| $k_{a,2}$ | association rate parameter for second binding | $RU^{-1}s^{-1}$ |
See above for discussion.
#### Differential Equation
$$\begin{align}
\frac{dR_1}{dt} &= 2k_{a,1}C(R_{max} - R_1 - 2R_2) - k_{d,1}R_1 - k_{a,2}R_1(R_{max} - R_1 - 2R_2) + 2k_{d,2}R_2 \\
\frac{dR_2}{dt} &= 2k_{a,2}R_1(R_{max} - R_1 - 2R_2) - 2k_{d,2}R_2 \\
R&=R_1+R_2
\end{align}$$
This is the equation used in the data simulation program.
#### Two Reference Versions
Here are two reference versions of this model.
##### Version 1: From 
Source: [[Huhn, Dushek 2024 - The molecular reach of antibodies]]
- [ ] Check publication date - seems to be 2025 #BS2026
Where:
$[Ab]$: Concentration of Analyte (antibody)
$A$: Concentration of unbound, immobilized receptor (antigen)
$B$: Concentration of antibody-antigen complex, single-bound
$C$: Concentration of antibody-antigent complex, double-bound

$$\begin{align}
\frac{dA}{dt} &= -2k_{on}[Ab]A + k_{off}B - k_{on,b}AB + 2k_{off}C \\
\frac{dB}{dt} &= 2k_{on}[Ab]A - k_{off}B - k_{on,b}AB + 2k_{off}C \\
\frac{dC}{dt} &= k_{on,b}AB - 2k_{off}C
\end{align}$$
This model is close to what we use in the simulation program. I standardized the variables and parameters, then eliminated one equation by using
$$\text{Available} = \text{Total} - B - 2C$$
##### Version 2: 
As recounted in Biacore Manuals:
$$\begin{align}
\frac{ds_1}{dt} &= k_{a1}A(s_{max} - s_1 - 2s_2) - k_{d1}s_1 - k_{a2}s_1(s_{max} - s_1 - 2s_2) + 2k_{d2}s_2 \\  
\frac{ds_2}{dt} &= k_{a2}s_1(s_{max} - s_1 - 2s_2) - 2k_{d2}s_2 \\
s &= s_1 + s_2
\end{align}$$
- This model eliminates a factor of two (absorbed into the parameters) and assumes a single value for $k_d$.
##### Notes
I retained the factor of two, and retained two separate values for $k_d$. 
- [ ] Justify this in future update #BS2026 
# Appendix / Notes / Scratch/Temporary

## To Do

Adjust Bivalent analyte model - one kD value only
Work on "Apparent KD" values derivations
