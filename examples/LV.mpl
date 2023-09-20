printf("\n================================================\n");
printf("= Software demo for the Lotka-Volterra example =\n");
printf("================================================\n\n");

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):


with(ListTools):
with(ArrayTools):
with(VectorCalculus):
with(LinearAlgebra):


#setting up the system
sys_rhs := [
            a*X[1](t) -b*X[1](t)*X[2](t),
            -c*X[2](t) + d*X[1](t)*X[2](t),
           y(t) - X[1](t)
    ]:   
params := [a, b, c, d]:
states := [seq(X[i], i=1..2)]:
outputs := [y]:
output_func := [X[1](t)]:
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

alg_subs := {seq(outputs_maxorders[j][2] = Y[j,0], j = 1 .. numelems(outputs)),
                 seq(seq(diff(outputs_maxorders[j][2], t$i) = Y[j, i], 
                         i = 1 .. degree(outputs_maxorders[j][1])
                    ), j = 1 .. numelems(outputs)
                 )
            }:
eq_alg := expand(subs(alg_subs, IOeqs)):

printf("\n============================================\n"):
printf("IO-identifiable functions of parameters are:\n\n"):
IO_coeffs := [seq(i, i in coeffs(expand(IOeqs[1]), map(lhs, alg_subs)))]:
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

num_out := 1:
printf("\n======================================\n");
printf("Computation for choosing alpha-tildes:\n");
for ov in LieDer do
   LieDerNo_t := subs([seq(xv(t) = xv, xv in states)], ov);
   Jac := Jacobian(LieDerNo_t, states):
   size_minor := degree(outputs_maxorders[num_out][1]):
   Q := Determinant(Jac[1..size_minor, 1..size_minor]):
   
   printf("\nA non-zero maximal minor of the Jacobian is\n\nQ = %a\n\n", Q);
   
   printf("Its numerator is\n\n%a\n\n", numer(Q)):
   
   Q0 := [coeffs(Q, states)][-1]:
   
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

chosen_tildes := [a_tilde = a, b_tilde = 1, c_tilde = c, d_tilde = d]:

for ch in chosen_tildes do printf("%a\n", ch); od;
wvars := [seq(w[j], j = 1..numelems(states))]:
x_to_wvars := [seq(states[i] = wvars[i], i = 1..numelems(states))]:
to_chosen_tildes := [
                      seq(params[i] 
                           = rhs(chosen_tildes[i]), 
                          i = 1..numelems(params))
                     ]:

printf("\n=============================\n");
printf("The resulting ODE systems is:\n\n");

for i from 1 to numelems(sys_rhs)-1 do 
   new_diffeq[i] := subs(to_chosen_tildes, simplify(sys_rhs[i], symbolic)): 
   printf("%a' = %a\n", wvars[i], subs(x_to_wvars, new_diffeq[i]));
od:
# Constructing the polynomial system for the change of variables
no_t_states := [seq(states[i](t) = states[i], i = 1..numelems(states))]:
LieDer_no_t := subs(no_t_states, LieDer):
LieDerW := subs([op(no_t_states), op(x_to_wvars), op(to_chosen_tildes)], LieDer_no_t):

for i from 1 to numelems(wvars) do 
   eq_conv[i] := LieDer_no_t[1][i] = LieDerW[1][i] 
od:
conv_sols := simplify(solve([seq(eq_conv[i],i=1..numelems(wvars))], wvars)):
printf("\n=========================================\n");
printf("The corresponding change of variables is:\n\n");

for s in conv_sols[1] do printf("%a\n", s): od;
