printf("\n=====================================\n");
printf("= Software demo for the CRN example =\n");
printf("=====================================\n\n");

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):


with(ListTools):
with(ArrayTools):
with(VectorCalculus):
with(LinearAlgebra):


#setting up the system
sys_rhs := [
          - X(t) * AUX(t) * k1 - X(t) * AXU(t) * k1 - 2 * X(t) * AUU(t) * k1 + AUX(t) * k2 + 2 * AXX(t) * k2 + AXU(t) * k2,
           -2 * X(t) * AUU(t) * k1 + AUX(t) * k2 + AXU(t) * k2,
           -X(t) * AUX(t) * k1 + X(t) * AUU(t) * k1 - AUX(t) * k2 + AXX(t) * k2,
            X(t) * AUX(t) * k1 + X(t) * AXU(t) * k1 - 2 * AXX(t) * k2,
           -X(t) * AXU(t) * k1 + X(t) * AUU(t) * k1 + AXX(t) * k2 - AXU(t) * k2
    ]:
params := [k1, k2]:
states := [X, AUU, AUX, AXX, AXU]:
outputs := [y]:
output_func := [X(t)]:
syst := [
         seq(diff(states[i](t),t) - sys_rhs[i], i = 1..numelems(states)), 
         seq(outputs[i](t) - output_func[i], i = 1..numelems(outputs))
        ]:
printf("==================\n");
printf("The ODE system is:\n\n");

no_t_states := [seq(states[i](t) = states[i], i = 1..numelems(states))]:
for i from 1 to numelems(states) do printf("%a' = %a\n", states[i],  subs(no_t_states, sys_rhs[i])) od;
for i from 1 to numelems(outputs) do printf("%a\n", outputs[i] = subs(no_t_states, output_func[i])) od;
# Computing input-output equations
R := DifferentialRing(
                      blocks = [states, outputs], 
                      derivations = [t], 
                      arbitrary = params
                     ):
eq := Equations(RosenfeldGroebner(syst, R))[1]:
IOeqs := simplify(eq[-numelems(outputs)..-1]):

printf("\n=====================\n");
printf("The IO-equations are:\n\n");

for eq in IOeqs do printf("%a\n\n", simplify(eq)); od;


# Replacing derivatives with regular variables for further computation

Rout := DifferentialRing(
                      blocks = [outputs], 
                      derivations = [t], 
                      arbitrary = params
                     ):
LD := LeadingDerivative(IOeqs, Rout):

outputs_maxorders := [seq([FactorDerivative(v, R)], v in LD)]:

alg_subs := {seq(outputs_maxorders[j][2] = Y[j,0], j = 1..numelems(outputs)),
                 seq(seq(diff(outputs_maxorders[j][2], t$i) = Y[j, i], 
                         i = 1 .. degree(outputs_maxorders[j][1])
                    ), j = 1 .. numelems(outputs)
                 )
            }:
eq_alg := expand(subs(alg_subs, IOeqs)):

printf("\n============================================\n"):
printf("IO-identifiable functions of parameters are:\n\n"):

IO_coeffs := [seq(i, i in coeffs(IOeqs[1], map(lhs, alg_subs)))]:
IO_coeffs_nonumbers := []:

for i in IO_coeffs do 
  if not type(i, integer) then IO_coeffs_nonumbers := [op(IO_coeffs_nonumbers), i] fi; 
od;

for i in IO_coeffs_nonumbers do 
  printf("%a\n", i) 
od;
# Compute Lie derivatives of the y-functions that participated in the IO-equations
R2 := DifferentialRing(
                      blocks = [outputs, states], 
                      derivations = [t], 
                      arbitrary = params
                     ):
eq2 := RosenfeldGroebner(syst, R2):

LieDer := [seq([op(NormalForm(outputs_maxorders[j][2], eq2)),
                 seq(op(NormalForm(diff(outputs_maxorders[j][2], t$i), eq2)), 
                         i = 1 .. degree(outputs_maxorders[j][1]))]
                     , j = 1 .. numelems(outputs)
                 )
            ]:

wvars := [seq(w[i], i = 1..3)]:

reduction_subs := [
                  X(t) = wvars[1], 
                  AUU(t) = wvars[2], 
                  AUX(t) = wvars[3], 
                  AXX(t) = 0, 
                  AXU(t) = 0
                  ]:

printf("\n=====================================\n");
printf("The chosen reduction substitution is:\n\n");

for sb in reduction_subs do 
  printf("%a\n\n", subs(no_t_states, simplify(sb))); 
od;

LieDerW := subs(reduction_subs, LieDer):
eqs_WX := [seq(LieDer[1][i] = LieDerW[1][i], i = 1..numelems(LieDer[1]))]:

ChangeVars := solve(eqs_WX, wvars)[1]:

printf("=========================================\n");
printf("The corresponding change of variables is:\n\n");
for sb in ChangeVars do 
  printf("%a\n\n", subs(no_t_states, simplify(sb))); 
od;
h := LieDerW:
alg_subs_rev := {seq(rhs(ex) = lhs(ex), ex in alg_subs)}:
# Computing the state-space form ODEs for the f's calculated above

# solving a linear system to find the new ODE system
dH := Concatenate(1, seq(Jacobian(h[i][1..-2], wvars), i = 1..numelems(outputs))):
H := Concatenate(1, seq(convert(h[i][2..-1], Vector), i = 1..numelems(outputs))):

reparam :=  convert(simplify(LinearSolve(dH, H)), list):

printf("==========================\n");
printf("The reduced ODE system is:\n\n");

for i from 1 to numelems(reparam) do 
   new_diffeq[i] := simplify(reparam[i], symbolic): 
   printf("%a' = %a\n",wvars[i], new_diffeq[i]);
od:
# Checking the result by recomputing the input-output equations
# and making sure that these are the same as the original input-output equations

w_t_subs := {seq(w[i] = w[i](t), i = 1..numelems(reparam))}:
for j from 1 to numelems(reparam) do 
   new_diffeq_t[j] := subs(w_t_subs, new_diffeq[j]) 
od:
ChangeVarsRev := [
                  seq(rhs(ChangeVars[i]) 
                         = lhs(ChangeVars[i]), 
                      i = 1..numelems(ChangeVars))
                 ]:
outputequations := subs(ChangeVarsRev, syst[-numelems(outputs)..-1]):
outputequations := subs(w_t_subs, outputequations):
syst_new := [seq(
                  diff(w[i](t), t) - new_diffeq_t[i], 
                  i = 1..numelems(reparam)
                ),
             op(outputequations)  
            ]:
R2 := DifferentialRing(
                      blocks = [wvars, outputs], 
                      derivations = [t], 
                      arbitrary = params
                     ):
eq := Equations(RosenfeldGroebner(syst_new, R2))[1]:
printf("\n==================================\n");
printf("The new IO-equations are the same:\n\n");
for i from 1 to numelems(outputs) do 
     printf("%a\n", evalb(simplify(eq[-i] - IOeqs[numelems(outputs)-i+1]=0))): 
od;
num_out := 1:
printf("\n======================================\n");
printf("Computation for choosing alpha-tildes:\n");
for ov in LieDerW do
   LieDerNo_t := subs([seq(xv(t) = xv, xv in wvars)], ov);
   Jac := Jacobian(LieDerNo_t, wvars):
   size_minor := degree(outputs_maxorders[num_out][1]):
   Q := Determinant(Jac[1..size_minor, 1..size_minor]):
   
   printf("\nA non-zero maximal minor of the Jacobian is\n\nQ = %a\n\n", Q);
   
   printf("Its numerator is\n\n%a\n\n", numer(Q)):
   
   Q0 := [coeffs(Q,wvars)][-1]:
   
   printf("Which has non-zero coefficient\n\nQ0 = %a\n", Q0);
   num_out += 1:
od:

params_tilde := [seq(cat(p,_tilde), p in params)]:
params_sub := [ 
               seq(params[i] 
                      = params_tilde[i], 
                   i = 1..numelems(params))
              ]:
tilde_eqs := [seq(i = subs(params_sub, i), i in IO_coeffs_nonumbers)]:
tilde_system := [op(tilde_eqs), subs(params_sub, Q0) <> 0]:

printf("\nThe corresponding system to choose the tilde-parameters is:\n\n");

for s in tilde_system do printf("%a\n", s) od;

printf("\nWe pick the following solution:\n\n");

chosen_tildes := [k1_tilde = k1, k2_tilde = 1]:

for c in chosen_tildes do printf("%a\n", c); od;
vvars := [seq(v[j], j = 1..numelems(wvars))]:
w_to_vvars := [seq(wvars[i] = vvars[i], i = 1..numelems(wvars))]:
to_chosen_tildes := [
                      seq(params[i] 
                           = rhs(chosen_tildes[i]), 
                          i = 1..numelems(params))
                     ]:

printf("\n=============================\n");
printf("The resulting ODE systems is:\n\n");

for i from 1 to numelems(reparam) do 
   new_diffeq[i] := subs(to_chosen_tildes, simplify(reparam[i], symbolic)): 
   printf("%a' = %a\n", vvars[i], subs(w_to_vvars, new_diffeq[i]));
od:
# Constructing the polynomial system for the change of variables
LieDerV := subs([op(w_to_vvars), op(to_chosen_tildes)], LieDerW):

for i from 1 to numelems(vvars) do 
   eq_conv[i] := LieDerW[1][i] = LieDerV[1][i] 
od:
conv_sols := simplify(solve([seq(eq_conv[i],i=1..numelems(vvars))], vvars)):

printf("\n=========================================\n");
printf("The corresponding change of variables is:\n\n");

for s in conv_sols[1] do printf("%a\n", s): od;

composed_change := subs(
                        [seq(states[i](t) = states[i], i=1..numelems(states))], 
                        simplify(subs(ChangeVars, conv_sols[1]))
                       ):

printf("\n==========================================\n");
printf("The final composed change of variables is:\n\n");
for s in composed_change do printf("%a\n", s): od;
