
restart:
printf("\n========================================\n");
printf("= Software demo for the bilinear model =\n");
printf("========================================\n\n");

with(PolynomialIdeals):
FieldToIdeal := proc(gens)
    # Input: generators of a subfield of the field of rational functions
    # Computes the MSQ ideal of the field with the new variables of the form x_aux
    # See: https://mediatum.ub.tum.de/doc/685465/document.pdf Definition 2.16
    #      https://doi.org/10.1016/j.jsc.2005.09.010          Lemma 2.1
    local all_vars, subs_dupl, all_dupl, common_denom, polys, f, gb:
    all_vars := indets(gens):
    subs_dupl := map(v -> v = cat(v, _aux), all_vars):
    all_dupl := sort([op(map(v -> subs(subs_dupl, v), all_vars))]):
    common_denom := 1:
    polys := []:
    for f in gens do
        common_denom := lcm(common_denom, denom(f)):
        polys := [op(polys), numer(f) * subs(subs_dupl, denom(f)) - subs(subs_dupl, numer(f)) * denom(f)]:
    end do:
    gb := Groebner[Basis]([op(polys), subs(subs_dupl, common_denom) * t - 1], tdeg(t, op(all_dupl))):
    gb := Groebner[Walk](gb, tdeg(t, op(all_dupl)), lexdeg([t], [op(all_dupl)])):
    gb := select(p -> not (t in indets(p)), gb):
    return PolynomialIdeal(gb, variables=all_dupl):
end proc:


#------------------------------------------------------------------------------

FilterGenerators := proc(ideal)
    # Input: ideal over a rational function field
    # Computes a simplified set of generators of the field of definition
    local gb, gens, p, cf, lc, gsorted, result, big_ideal, cur_ideal, g, new_ideal:
    gb := Groebner[Basis](ideal, tdeg(op(IdealInfo[Variables](ideal)))):
    gens := {}:
    for p in gb do
        cf := sort([coeffs(p, IdealInfo[Variables](ideal))], (a, b) -> length(convert(a, string)) < length(convert(b, string))):
        lc := cf[1]:
        cf := map(c -> c / lc, cf):
        gens := {op(gens), op(cf[2..nops(cf)])}:
    end do:
    gsorted := sort([op(gens)], (a, b) -> length(convert(a, string)) < length(convert(b, string)));
    result := {}:
    big_ideal := FieldToIdeal(gens):
    cur_ideal := FieldToIdeal(result):
    for g in gsorted do
        if big_ideal = cur_ideal then
            return result:
        end if:
        new_ideal := FieldToIdeal({op(result), g}):
        if new_ideal <> cur_ideal then
            # a dirty hack to transform -1/a to a
            if convert(g, string)[1] = "-" then
                g := -g:
            end:
            if convert(g, string)[1] = "1" then
                g := 1 / g:
            end:
            result := {op(result), g}:
            cur_ideal := new_ideal:
        end if:
    end do:
    return result:
end proc:

# Loading packages
with(DifferentialAlgebra):
with(DifferentialAlgebra[Tools]):


with(ListTools):
with(ArrayTools):
with(VectorCalculus):
with(LinearAlgebra):


#setting up the system
sys_rhs := [
            - p1 * x1(t) + p2 * u(t),
            - p3 * x2(t) + p4 * u(t),
            - (p1 + p3) * x3(t) + (p4 * x1(t) + p2 * x2(t)) * u(t) 
    ]:   
params := [p1, p2, p3, p4]:
states := [x1, x2, x3]:
outputs := [y]:
output_func := [x3(t)]:
inputs := [u]:
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
IO_coeffs := [seq(i, i in coeffs(expand(IOeqs[1]), map(lhs, alg_subs)))]:
IO_coeffs := FilterGenerators(FieldToIdeal(IO_coeffs)):
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

num_out := 1:
printf("\n======================================\n");
printf("Computation for choosing alpha-tildes:\n");
for ov in LieDer do
   LieDerNo_t := subs([seq(xv(t) = xv, xv in [op(states), op(inputs)])], ov);
   Jac := simplify(Jacobian(LieDerNo_t, states));
   size_minor := degree(outputs_maxorders[num_out][1]):
   Q := Determinant(Jac[1..size_minor, 1..size_minor]):
   
   printf("\nA non-zero maximal minor of the Jacobian is\n\nQ = %a\n\n", Q);
   
   printf("Its numerator is\n\n%a\n\n", numer(Q)):
   
   Q0 := [coeffs(expand(numer(Q)), [op(states), op(inputs)])][-1]:
   
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

chosen_tildes := [                   p1_tilde = p1,                    p2_tilde = 1,                    p3_tilde = p3,                    p4_tilde = p2 * p4                 ]:

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
LieDer_no_t := subs([op(no_t_states), op(no_t_inputs)], LieDer):
LieDerW := subs([op(no_t_states), op(x_to_wvars), op(to_chosen_tildes)], LieDer_no_t):

for i from 1 to numelems(wvars) do 
   eq_conv[i] := LieDer_no_t[1][i] = LieDerW[1][i] 
od:
conv_sols := simplify(solve([seq(eq_conv[i],i=1..numelems(wvars))], wvars)):
printf("\n=========================================\n");
printf("The corresponding change of variables is:\n\n");

for s in conv_sols[1] do printf("%a\n", s): od;

