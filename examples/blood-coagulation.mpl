
restart:
printf("\n================================================\n");
printf("= Software demo for the blood coagulation model =\n");
printf("================================================\n\n");

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):
with(ListTools):
with(ArrayTools):
with(VectorCalculus):
with(LinearAlgebra):
read("common.mpl"):

#setting up the system
sys_rhs := [
                k1 * beta - h1 * IXa(t),
                k2 * IIa(t) + k3 * APC(t) * VIIIa(t)/(b1 + VIIIa(t)) - h2 * VIIIa(t),
                k5 * IXa(t) * VIIIa(t)/(b2 + VIIIa(t)) - h3 * Xa(t),
                k6 * IIa(t) - k7 * APC(t) * Va(t)/(b3 + Va(t)) - h4 * Va(t),
                k8 * IIa(t) - h5 * APC(t),
                k9 * Xa(t) * Va(t)/(b4 + Va(t)) - h6 * IIa(t)
    ]:   
params := [           b1, b2, b3, b4, beta,            h1, h2, h3, h4, h5, h6,            k1, k2, k3, k5, k6, k7, k8, k9          ]:
states := [IXa, VIIIa, Xa, Va, APC, IIa]:
outputs := [y1, y2, y3, y4]:
output_func := [Xa(t), IIa(t), IXa(t), APC(t)]:
inputs := []:
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
                      blocks = [states, outputs, inputs], 
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
                      blocks = [outputs,inputs], 
                      derivations = [t], 
                      arbitrary = params
                     ):
LD := LeadingDerivative(IOeqs, Rout):
outputs_maxorders := [seq([FactorDerivative(v, R)], v in LD)]:
Rin := DifferentialRing(
                      blocks = [inputs,outputs], 
                      derivations = [t], 
                      arbitrary = params
                     ):
LD := LeadingDerivative(IOeqs, Rin):
inputs_maxorders := [seq([FactorDerivative(v, R)], v in LD)]:

alg_subs := {seq(outputs_maxorders[j][2] = Y[j,0], j = 1..numelems(outputs)),
                 seq(seq(diff(outputs_maxorders[j][2], t$i) = Y[j,i], 
                         i = 1 .. degree(outputs_maxorders[j][1])
                    ), j = 1 .. numelems(outputs)
                 ),
              seq(inputs_maxorders[j][2] = U[j,0], j = 1..numelems(inputs)),
                 seq(seq(diff(inputs_maxorders[j][2], t$i) = U[j,i], 
                         i = 1 .. degree(inputs_maxorders[j][1])
                    ), j = 1 .. numelems(inputs)
                 )
            }:
eq_alg := expand(subs(alg_subs, IOeqs)):

printf("\n============================================\n"):
printf("IO-identifiable functions of parameters are:\n\n"):
IO_coeffs := [seq(i, i in coeffs(expand(IOeqs[1]/(b1 * k2)), map(lhs, alg_subs)))]:
IO_coeffs := [op(IO_coeffs), seq(i, i in coeffs(expand(IOeqs[2]/(b3 * k6)), map(lhs, alg_subs)))]:
IO_coeffs := [op(IO_coeffs), seq(i, i in coeffs(expand(IOeqs[3]), map(lhs, alg_subs)))]:
IO_coeffs := [op(IO_coeffs), seq(i, i in coeffs(expand(IOeqs[4]), map(lhs, alg_subs)))]:
IO_coeffs := SimplifyRationalFunctions(IO_coeffs):
IO_coeffs_nonumbers := []:

for i in IO_coeffs do 
  if not type(i, integer) then IO_coeffs_nonumbers := [op(IO_coeffs_nonumbers), i] fi; 
od;

for i in IO_coeffs_nonumbers do 
  printf("%a\n", i) 
od;
# Compute Lie derivatives of the y-functions that participated in the IO-equations
R2 := DifferentialRing(
                      blocks = [outputs, states, inputs], 
                      derivations = [t], 
                      arbitrary = params
                     ):
eq2 := RosenfeldGroebner(syst, R2):

LieDer := simplify([seq([op(NormalForm(outputs_maxorders[j][2], eq2)),
                 seq(op(NormalForm(diff(outputs_maxorders[j][2], t$i), eq2)), 
                         i = 1 .. degree(outputs_maxorders[j][1]))]
                     , j = 1 .. numelems(outputs)
                 )
            ]):


size_minor := add(degree(outputs_maxorders[num_out][1]), num_out = 1..numelems(outputs)):

printf("\n======================================\n");
printf("Computation for choosing alpha-tildes:\n");
LieDerNo_t := subs([seq(xv(t) = xv, xv in [op(states), op(inputs)])], LieDer):
Jac := simplify(Jacobian(Flatten(LieDerNo_t), states)):
  
Q := Determinant(Jac[1..size_minor, 1..size_minor]):
   
printf("\nA non-zero maximal minor of the Jacobian is\n\nQ = %a\n\n", Q);
   
printf("Its numerator is\n\n%a\n\n", numer(Q)):
   
Q0 := [coeffs(expand(numer(Q)), [op(states), op(inputs)])][-1]:
   
printf("Which has non-zero coefficient\n\nQ0 = %a\n", Q0);

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

chosen_tildes := [                  b1_tilde = b1 * k3_tilde / k3,                   b2_tilde = b2 * k3_tilde / k3,                   b3_tilde = b3 * k7_tilde / k7,                   b4_tilde = b4 * k7_tilde / k7,                   beta_tilde = k1 * beta / k1_tilde,                   h1_tilde = h1,                   h2_tilde = h2,                   h3_tilde = h3,                   h4_tilde = h4,                   h5_tilde = h5,                   h6_tilde = h6,                   k1_tilde = k1_tilde,                   k2_tilde = k3_tilde * k2 / k3,                   k3_tilde = k3_tilde,                   k5_tilde = k5,                   k6_tilde = k7_tilde * k6 / k7,                   k7_tilde = k7_tilde,                   k8_tilde = k8,                   k9_tilde = k9                 ]:
# updating the chosen solution with specific values for the free variables
chosen_tildes := subs([k1_tilde = 1, k3_tilde = 1, k7_tilde = 1], chosen_tildes):
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

for i from 1 to numelems(sys_rhs) do 
   new_diffeq[i] := subs(to_chosen_tildes, simplify(sys_rhs[i], symbolic)): 
   printf("%a' = %a\n", wvars[i], subs(x_to_wvars, new_diffeq[i]));
od:
# Constructing the polynomial system for the change of variables
no_t_states := [seq(states[i](t) = states[i], i = 1..numelems(states))]:
no_t_inputs := [seq(inputs[i](t) = inputs[i], i = 1..numelems(inputs))]:
LieDer_no_t := Flatten(subs([op(no_t_states), op(no_t_inputs)], LieDer)):
LieDerW := Flatten(subs([op(no_t_states), op(x_to_wvars), op(to_chosen_tildes)], LieDer_no_t)):

for i from 1 to numelems(wvars)+2 do 
   eq_conv[i] := LieDer_no_t[i] = LieDerW[i] 
od:
conv_sols := simplify(solve([seq(eq_conv[i],i=1..numelems(wvars))], wvars)):
printf("\n=========================================\n");
printf("The corresponding change of variables is:\n\n");

for s in conv_sols[1] do printf("%a\n", s): od;
