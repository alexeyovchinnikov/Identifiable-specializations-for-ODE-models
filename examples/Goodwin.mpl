
restart:
printf("\n====================================================\n");
printf("= Software demo for the Goodwin Oscillator example =\n");
printf("====================================================\n\n");

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
            -b * x1(t) + 1 / (c + x4(t)),
            alpha * x1(t) - beta * x2(t),
            gama * x2(t) - delta * x3(t),
            sigma * x4(t) * (gama * x2(t) - delta * x3(t)) / x3(t)
    ]:   
params := [b, c, alpha, beta, gama, delta, sigma]:
states := [x1, x2, x3, x4]:
outputs := [y1]:
output_func := [x1(t)]:
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
                      blocks = [outputs, states], 
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

num_out := 1:
printf("\n======================================\n");
printf("Computation for choosing alpha-tildes:\n");
for ov in LieDer do
   LieDerNo_t := subs([seq(xv(t) = xv, xv in states)], ov);
   Jac := simplify(Jacobian(LieDerNo_t, states));
   size_minor := degree(outputs_maxorders[num_out][1]):
   
   printf("\nThe minor has size:\n"):
   printf("%a\n", size_minor):
   Q := Determinant(Jac[[1..size_minor], [1..size_minor]]):
   
   printf("\nA non-zero maximal minor of the Jacobian is\n\nQ = %a\n\n", Q);
   
   printf("Its numerator is\n\n%a\n\n", numer(Q)):
   
   Q0 := [coeffs(expand(numer(Q)), states)][-1]:
   
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

chosen_tildes := [                  b_tilde = b,                   c_tilde = c,                   alpha_tilde = 1 ,                   beta_tilde = beta,                   gama_tilde = 1,                   delta_tilde = delta,                   sigma_tilde = sigma                  ]:

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
   printf("%a' = %a\n", wvars[i], subs(x_to_wvars, simplify(new_diffeq[i])));
od:
# Constructing the polynomial system for the change of variables
no_t_states := [seq(states[i](t) = states[i], i = 1..numelems(states))]:
LieDer_no_t := subs(no_t_states, LieDer):
LieDerW := subs([op(no_t_states), op(x_to_wvars), op(to_chosen_tildes)], LieDer_no_t):

for i from 1 to numelems(LieDer_no_t[1]) do 
   eq_conv[i] := LieDer_no_t[1][i] = LieDerW[1][i] 
od:

conv_sols := simplify(solve([seq(eq_conv[i],i=1..numelems(LieDer_no_t[1]))], wvars)):
printf("\n=========================================\n");
printf("The corresponding change of variables is:\n\n");

for s in conv_sols[1] do printf("%a\n", simplify(s)): od;
