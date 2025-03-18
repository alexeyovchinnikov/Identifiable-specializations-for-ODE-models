

restart:
printf("\n============================================\n");
printf("= Software demo for the Akt - Fujita model =\n");
printf("============================================\n\n");

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
            
               X1_t * u(t) + X9(t) * r_1_k2 -
               X1(t) * X1_t - X9(t) * r_1_k1,
         
               X9(t) * r_9_k1 - X2(t) * r_4_k1 +
               X3(t) * r_2_k2 +
               X3(t) * r_3_k1 - X4(t) * X2(t) * r_2_k1,
          
               X4(t) * X2(t) * r_2_k1 - X3(t) * r_3_k1 -
               X3(t) * r_2_k2,
         
               X5(t) * r_7_k1 + X3(t) * r_2_k2 -
               X4(t) * X2(t) * r_2_k1,
          
               X7(t) * r_5_k2 - X5(t) * r_7_k1 +
               X7(t) * r_6_k1 +
               X3(t) * r_3_k1 - X6(t) * X5(t) * r_5_k1,
           
               X7(t) * r_5_k2 + X8(t) * r_8_k1 -
               X6(t) * X5(t) * r_5_k1,
           
               X6(t) * X5(t) * r_5_k1 - X7(t) * r_6_k1 -
               X7(t) * r_5_k2,
            X7(t) * r_6_k1 - X8(t) * r_8_k1,
           
               X9(t) * r_1_k1 - X9(t) * r_9_k1 -
               X9(t) * r_1_k2
    ]:   
params := [X1_t, a1, a2, a3, r_1_k1, r_1_k2, r_2_k1, r_2_k2, r_3_k1, r_4_k1, r_5_k1, r_5_k2, r_6_k1, r_7_k1, r_8_k1, r_9_k1]:
states := [X1, X2, X3, X4, X5, X6, X7, X8, X9]:
outputs := [y1, y2, y3]:
output_func := [a1 * (X2(t) + X3(t)), a2 * (X5(t) + X7(t)), a3 * X8(t)]:
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

# The IO-iquations below were calculated using StructuralIdentifiability.jl in Julia using the function "find_ioequations"   

eqs := [a2^2*a3*r_6_k1*r_8_k1^2*y3(t)*diff(y3(t), t$2) - a2^2*a3*r_6_k1*r_8_k1^2*diff(y3(t), t)^2 + a2^2*a3*r_6_k1*r_8_k1*y3(t)*diff(y3(t), t$3) - a2^2*a3*r_6_k1*r_8_k1*diff(y3(t), t)*diff(y3(t), t$2) + a2^2*a3*r_6_k1*diff(y3(t), t)*diff(y3(t), t$3) - a2^2*a3*r_6_k1*diff(y3(t), t$2)^2 - a2^2*r_5_k1*r_6_k1*r_8_k1^2*y3(t)^2*diff(y3(t), t) - 2*a2^2*r_5_k1*r_6_k1*r_8_k1*y3(t)*diff(y3(t), t)^2 - a2^2*r_5_k1*r_6_k1*diff(y3(t), t)^3 - a2^2*r_5_k1*r_8_k1^3*y3(t)^2*diff(y3(t), t) - a2^2*r_5_k1*r_8_k1^2*y3(t)^2*diff(y3(t), t$2) - 2*a2^2*r_5_k1*r_8_k1^2*y3(t)*diff(y3(t), t)^2 - 2*a2^2*r_5_k1*r_8_k1*y3(t)*diff(y3(t), t)*diff(y3(t), t$2) - a2^2*r_5_k1*r_8_k1*diff(y3(t), t)^3 - a2^2*r_5_k1*diff(y3(t), t)^2*diff(y3(t), t$2) - a2*a3^2*r_5_k2*r_6_k1^2*r_8_k1*y2(t)*diff(y3(t), t) + a2*a3^2*r_5_k2*r_6_k1^2*r_8_k1*y3(t)*diff(y2(t), t) - a2*a3^2*r_5_k2*r_6_k1^2*y2(t)*diff(y3(t), t$2) + a2*a3^2*r_5_k2*r_6_k1^2*diff(y2(t), t)*diff(y3(t), t) - a2*a3^2*r_6_k1^3*r_8_k1*y2(t)*diff(y3(t), t) + a2*a3^2*r_6_k1^3*r_8_k1*y3(t)*diff(y2(t), t) - a2*a3^2*r_6_k1^3*y2(t)*diff(y3(t), t$2) + a2*a3^2*r_6_k1^3*diff(y2(t), t)*diff(y3(t), t) - a2*a3^2*r_6_k1^2*r_8_k1*y2(t)*diff(y3(t), t$2) + a2*a3^2*r_6_k1^2*r_8_k1*diff(y2(t), t)*diff(y3(t), t) - a2*a3^2*r_6_k1^2*y2(t)*diff(y3(t), t$3) + a2*a3^2*r_6_k1^2*diff(y2(t), t)*diff(y3(t), t$2) + 2*a2*a3*r_5_k1*r_6_k1^2*r_8_k1*y2(t)*y3(t)*diff(y3(t), t) + 2*a2*a3*r_5_k1*r_6_k1^2*y2(t)*diff(y3(t), t)^2 + 2*a2*a3*r_5_k1*r_6_k1*r_8_k1^2*y2(t)*y3(t)*diff(y3(t), t) + 2*a2*a3*r_5_k1*r_6_k1*r_8_k1*y2(t)*y3(t)*diff(y3(t), t$2) + 2*a2*a3*r_5_k1*r_6_k1*r_8_k1*y2(t)*diff(y3(t), t)^2 + 2*a2*a3*r_5_k1*r_6_k1*y2(t)*diff(y3(t), t)*diff(y3(t), t$2) - a3^2*r_5_k1*r_6_k1^3*y2(t)^2*diff(y3(t), t) - a3^2*r_5_k1*r_6_k1^2*r_8_k1*y2(t)^2*diff(y3(t), t) - a3^2*r_5_k1*r_6_k1^2*y2(t)^2*diff(y3(t), t$2),

a1*a2*r_1_k1*r_4_k1*r_7_k1*r_8_k1*y3(t) + a1*a2*r_1_k1*r_4_k1*r_7_k1*diff(y3(t), t) - a1*a2*r_1_k2*r_4_k1*r_7_k1*r_8_k1*y3(t) - a1*a2*r_1_k2*r_4_k1*r_7_k1*diff(y3(t), t) - a1*a2*r_4_k1*r_7_k1*r_8_k1*r_9_k1*y3(t) - a1*a2*r_4_k1*r_7_k1*r_8_k1*diff(y3(t), t) - a1*a2*r_4_k1*r_7_k1*r_9_k1*diff(y3(t), t) - a1*a2*r_4_k1*r_7_k1*diff(y3(t), t$2) - a1*a3*r_1_k1*r_4_k1*r_6_k1*r_7_k1*y2(t) - a1*a3*r_1_k1*r_4_k1*r_6_k1*diff(y2(t), t) + a1*a3*r_1_k2*r_4_k1*r_6_k1*r_7_k1*y2(t) + a1*a3*r_1_k2*r_4_k1*r_6_k1*diff(y2(t), t) + a1*a3*r_4_k1*r_6_k1*r_7_k1*r_9_k1*y2(t) + a1*a3*r_4_k1*r_6_k1*r_7_k1*diff(y2(t), t) + a1*a3*r_4_k1*r_6_k1*r_9_k1*diff(y2(t), t) + a1*a3*r_4_k1*r_6_k1*diff(y2(t), t$2) + a2*a3*r_1_k1*r_3_k1*r_4_k1*r_6_k1*y1(t) + a2*a3*r_1_k1*r_3_k1*r_6_k1*diff(y1(t), t) - a2*a3*r_1_k2*r_3_k1*r_4_k1*r_6_k1*y1(t) - a2*a3*r_1_k2*r_3_k1*r_6_k1*diff(y1(t), t) - a2*a3*r_3_k1*r_4_k1*r_6_k1*r_9_k1*y1(t) - a2*a3*r_3_k1*r_4_k1*r_6_k1*diff(y1(t), t) - a2*a3*r_3_k1*r_6_k1*r_9_k1*diff(y1(t), t) - a2*a3*r_3_k1*r_6_k1*diff(y1(t), t$2),


a1^2*a2^4*r_2_k1*r_7_k1^3*r_8_k1^4*y3(t)^3*diff(y3(t), t) + a1^2*a2^4*r_2_k1*r_7_k1^3*r_8_k1^3*y3(t)^3*diff(y3(t), t$2) + 3*a1^2*a2^4*r_2_k1*r_7_k1^3*r_8_k1^3*y3(t)^2*diff(y3(t), t)^2 + 3*a1^2*a2^4*r_2_k1*r_7_k1^3*r_8_k1^2*y3(t)^2*diff(y3(t), t)*diff(y3(t), t$2) + 3*a1^2*a2^4*r_2_k1*r_7_k1^3*r_8_k1^2*y3(t)*diff(y3(t), t)^3 + 3*a1^2*a2^4*r_2_k1*r_7_k1^3*r_8_k1*y3(t)*diff(y3(t), t)^2*diff(y3(t), t$2) + a1^2*a2^4*r_2_k1*r_7_k1^3*r_8_k1*diff(y3(t), t)^4 + a1^2*a2^4*r_2_k1*r_7_k1^3*diff(y3(t), t)^3*diff(y3(t), t$2) + a1^2*a2^4*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*r_8_k1^3*y3(t)^3*diff(y3(t), t) + 3*a1^2*a2^4*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y3(t), t)^2 + 3*a1^2*a2^4*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*r_8_k1*y3(t)*diff(y3(t), t)^3 + a1^2*a2^4*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*diff(y3(t), t)^4 + a1^2*a2^4*r_3_k1*r_5_k1*r_7_k1^2*r_8_k1^4*y3(t)^3*diff(y3(t), t) + a1^2*a2^4*r_3_k1*r_5_k1*r_7_k1^2*r_8_k1^3*y3(t)^3*diff(y3(t), t$2) + 3*a1^2*a2^4*r_3_k1*r_5_k1*r_7_k1^2*r_8_k1^3*y3(t)^2*diff(y3(t), t)^2 + 3*a1^2*a2^4*r_3_k1*r_5_k1*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y3(t), t)*diff(y3(t), t$2) + 3*a1^2*a2^4*r_3_k1*r_5_k1*r_7_k1^2*r_8_k1^2*y3(t)*diff(y3(t), t)^3 + 3*a1^2*a2^4*r_3_k1*r_5_k1*r_7_k1^2*r_8_k1*y3(t)*diff(y3(t), t)^2*diff(y3(t), t$2) + a1^2*a2^4*r_3_k1*r_5_k1*r_7_k1^2*r_8_k1*diff(y3(t), t)^4 + a1^2*a2^4*r_3_k1*r_5_k1*r_7_k1^2*diff(y3(t), t)^3*diff(y3(t), t$2) + a1^2*a2^3*a3^2*r_3_k1*r_5_k2*r_6_k1^2*r_7_k1^2*r_8_k1^2*y2(t)*y3(t)*diff(y3(t), t) - a1^2*a2^3*a3^2*r_3_k1*r_5_k2*r_6_k1^2*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y2(t), t) + a1^2*a2^3*a3^2*r_3_k1*r_5_k2*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y3(t), t$2) + a1^2*a2^3*a3^2*r_3_k1*r_5_k2*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)*diff(y3(t), t)^2 - 2*a1^2*a2^3*a3^2*r_3_k1*r_5_k2*r_6_k1^2*r_7_k1^2*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t) + a1^2*a2^3*a3^2*r_3_k1*r_5_k2*r_6_k1^2*r_7_k1^2*y2(t)*diff(y3(t), t)*diff(y3(t), t$2) - a1^2*a2^3*a3^2*r_3_k1*r_5_k2*r_6_k1^2*r_7_k1^2*diff(y2(t), t)*diff(y3(t), t)^2 + a1^2*a2^3*a3^2*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1^2*y2(t)*y3(t)*diff(y3(t), t) - a1^2*a2^3*a3^2*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y2(t), t) + a1^2*a2^3*a3^2*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y3(t), t$2) + a1^2*a2^3*a3^2*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)*diff(y3(t), t)^2 - 2*a1^2*a2^3*a3^2*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t) + a1^2*a2^3*a3^2*r_3_k1*r_6_k1^3*r_7_k1^2*y2(t)*diff(y3(t), t)*diff(y3(t), t$2) - a1^2*a2^3*a3^2*r_3_k1*r_6_k1^3*r_7_k1^2*diff(y2(t), t)*diff(y3(t), t)^2 - a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y2(t), t$2) + a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y3(t)*diff(y2(t), t)*diff(y3(t), t) + a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t$2) - 2*a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y3(t)*diff(y3(t), t)*diff(y2(t), t$2) + a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1*diff(y2(t), t)*diff(y3(t), t)^2 + a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1^2*diff(y2(t), t)*diff(y3(t), t)*diff(y3(t), t$2) - a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1^2*diff(y3(t), t)^2*diff(y2(t), t$2) - a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)^2*diff(y2(t), t$3) + 2*a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)*diff(y3(t), t)*diff(y2(t), t$2) - a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1^2*diff(y2(t), t)*diff(y3(t), t)^2 - 2*a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y3(t), t)*diff(y2(t), t$3) + 2*a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t$2)*diff(y3(t), t$2) - 2*a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1*diff(y2(t), t)*diff(y3(t), t)*diff(y3(t), t$2) + 2*a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1*diff(y3(t), t)^2*diff(y2(t), t$2) - a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*diff(y2(t), t)*diff(y3(t), t$2)^2 - a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*diff(y3(t), t)^2*diff(y2(t), t$3) + 2*a1^2*a2^3*a3^2*r_3_k1*r_6_k1^2*r_7_k1*diff(y3(t), t)*diff(y2(t), t$2)*diff(y3(t), t$2) - a1^2*a2^3*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*r_8_k1^3*y3(t)^3*diff(y2(t), t) - 3*a1^2*a2^3*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y2(t), t)*diff(y3(t), t) - 3*a1^2*a2^3*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t)^2 - a1^2*a2^3*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*diff(y2(t), t)*diff(y3(t), t)^3 - 3*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*r_8_k1^3*y2(t)*y3(t)^2*diff(y3(t), t) - a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*r_8_k1^3*y3(t)^3*diff(y2(t), t) - 3*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*r_8_k1^2*y2(t)*y3(t)^2*diff(y3(t), t$2) - 6*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*r_8_k1^2*y2(t)*y3(t)*diff(y3(t), t)^2 - 3*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*r_8_k1^2*y3(t)^2*diff(y2(t), t)*diff(y3(t), t) - 6*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*r_8_k1*y2(t)*y3(t)*diff(y3(t), t)*diff(y3(t), t$2) - 3*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*r_8_k1*y2(t)*diff(y3(t), t)^3 - 3*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t)^2 - 3*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*y2(t)*diff(y3(t), t)^2*diff(y3(t), t$2) - a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^3*diff(y2(t), t)*diff(y3(t), t)^3 - a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*r_8_k1^3*y3(t)^3*diff(y2(t), t$2) - 2*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*r_8_k1^3*y3(t)^2*diff(y2(t), t)*diff(y3(t), t) - 2*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y2(t), t)*diff(y3(t), t$2) - 3*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y3(t), t)*diff(y2(t), t$2) - 4*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y3(t)*diff(y2(t), t)*diff(y3(t), t)^2 - 4*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t)*diff(y3(t), t$2) - 3*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*r_8_k1*y3(t)*diff(y3(t), t)^2*diff(y2(t), t$2) - 2*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*r_8_k1*diff(y2(t), t)*diff(y3(t), t)^3 - 2*a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*diff(y2(t), t)*diff(y3(t), t)^2*diff(y3(t), t$2) - a1^2*a2^3*a3*r_2_k1*r_6_k1*r_7_k1^2*diff(y3(t), t)^3*diff(y2(t), t$2) - 3*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y2(t)*y3(t)^2*diff(y3(t), t) - 6*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y3(t), t)^2 - 3*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1^2*y2(t)*diff(y3(t), t)^3 - a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)^2*diff(y2(t), t)*diff(y3(t), t) - 2*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t)^2 - a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1*diff(y2(t), t)*diff(y3(t), t)^3 - 3*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*r_8_k1^3*y2(t)*y3(t)^2*diff(y3(t), t) - 3*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y2(t)*y3(t)^2*diff(y3(t), t$2) - 6*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y2(t)*y3(t)*diff(y3(t), t)^2 - 6*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y3(t), t)*diff(y3(t), t$2) - 3*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*r_8_k1*y2(t)*diff(y3(t), t)^3 - 3*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1^2*y2(t)*diff(y3(t), t)^2*diff(y3(t), t$2) - a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1*r_8_k1^3*y3(t)^2*diff(y2(t), t)*diff(y3(t), t) - a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1*r_8_k1^2*y3(t)^2*diff(y2(t), t)*diff(y3(t), t$2) - 2*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1*r_8_k1^2*y3(t)*diff(y2(t), t)*diff(y3(t), t)^2 - 2*a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t)*diff(y3(t), t$2) - a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1*r_8_k1*diff(y2(t), t)*diff(y3(t), t)^3 - a1^2*a2^3*a3*r_3_k1*r_5_k1*r_6_k1*r_7_k1*diff(y2(t), t)*diff(y3(t), t)^2*diff(y3(t), t$2) - a1^2*a2^2*a3^3*r_3_k1*r_5_k2*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)^2*diff(y3(t), t) + a1^2*a2^2*a3^3*r_3_k1*r_5_k2*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y2(t), t) - a1^2*a2^2*a3^3*r_3_k1*r_5_k2*r_6_k1^3*r_7_k1^2*y2(t)^2*diff(y3(t), t$2) + a1^2*a2^2*a3^3*r_3_k1*r_5_k2*r_6_k1^3*r_7_k1^2*y2(t)*diff(y2(t), t)*diff(y3(t), t) - a1^2*a2^2*a3^3*r_3_k1*r_5_k2*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t) + a1^2*a2^2*a3^3*r_3_k1*r_5_k2*r_6_k1^3*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)^2 - a1^2*a2^2*a3^3*r_3_k1*r_5_k2*r_6_k1^3*r_7_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t$2) + a1^2*a2^2*a3^3*r_3_k1*r_5_k2*r_6_k1^3*r_7_k1*diff(y2(t), t)^2*diff(y3(t), t) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^4*r_7_k1^2*r_8_k1*y2(t)^2*diff(y3(t), t) + a1^2*a2^2*a3^3*r_3_k1*r_6_k1^4*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y2(t), t) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^4*r_7_k1^2*y2(t)^2*diff(y3(t), t$2) + a1^2*a2^2*a3^3*r_3_k1*r_6_k1^4*r_7_k1^2*y2(t)*diff(y2(t), t)*diff(y3(t), t) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^4*r_7_k1*r_8_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t) + a1^2*a2^2*a3^3*r_3_k1*r_6_k1^4*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)^2 - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^4*r_7_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t$2) + a1^2*a2^2*a3^3*r_3_k1*r_6_k1^4*r_7_k1*diff(y2(t), t)^2*diff(y3(t), t) + 2*a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y2(t), t$2) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y3(t)*diff(y2(t), t)^2 - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1^2*y2(t)*diff(y2(t), t)*diff(y3(t), t$2) + 2*a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1^2*y2(t)*diff(y3(t), t)*diff(y2(t), t$2) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1^2*diff(y2(t), t)^2*diff(y3(t), t) + 2*a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*y3(t)*diff(y2(t), t$3) - 2*a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*diff(y3(t), t)*diff(y2(t), t$2) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)*diff(y2(t), t$2) + a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1*r_8_k1*diff(y2(t), t)^2*diff(y3(t), t) + 2*a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1*y2(t)*diff(y3(t), t)*diff(y2(t), t$3) - 2*a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1*y2(t)*diff(y2(t), t$2)*diff(y3(t), t$2) + a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1*diff(y2(t), t)^2*diff(y3(t), t$2) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_7_k1*diff(y2(t), t)*diff(y3(t), t)*diff(y2(t), t$2) + a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_8_k1*y3(t)*diff(y2(t), t)*diff(y2(t), t$3) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*r_8_k1*y3(t)*diff(y2(t), t$2)^2 + a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*diff(y2(t), t)*diff(y3(t), t)*diff(y2(t), t$3) - a1^2*a2^2*a3^3*r_3_k1*r_6_k1^3*diff(y3(t), t)*diff(y2(t), t$2)^2 + 3*a1^2*a2^2*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y2(t)*y3(t)^2*diff(y2(t), t) + 6*a1^2*a2^2*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t) + 3*a1^2*a2^2*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*y2(t)*diff(y2(t), t)*diff(y3(t), t)^2 + 2*a1^2*a2^2*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)^2*diff(y2(t), t)^2 + 4*a1^2*a2^2*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)^2*diff(y3(t), t) + 2*a1^2*a2^2*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*diff(y2(t), t)^2*diff(y3(t), t)^2 + 3*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^3*r_8_k1^2*y2(t)^2*y3(t)*diff(y3(t), t) + 3*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^3*r_8_k1^2*y2(t)*y3(t)^2*diff(y2(t), t) + 3*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^3*r_8_k1*y2(t)^2*y3(t)*diff(y3(t), t$2) + 3*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^3*r_8_k1*y2(t)^2*diff(y3(t), t)^2 + 6*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^3*r_8_k1*y2(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t) + 3*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^3*y2(t)^2*diff(y3(t), t)*diff(y3(t), t$2) + 3*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^3*y2(t)*diff(y2(t), t)*diff(y3(t), t)^2 + 3*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y2(t)*y3(t)^2*diff(y2(t), t$2) + 4*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y2(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t) + 2*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y3(t)^2*diff(y2(t), t)^2 + 4*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t$2) + 6*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y3(t), t)*diff(y2(t), t$2) + 4*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t)^2 + 4*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y3(t)*diff(y2(t), t)^2*diff(y3(t), t) + 4*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*y2(t)*diff(y2(t), t)*diff(y3(t), t)*diff(y3(t), t$2) + 3*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*y2(t)*diff(y3(t), t)^2*diff(y2(t), t$2) + 2*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1^2*diff(y2(t), t)^2*diff(y3(t), t)^2 + 2*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)^2*diff(y2(t), t)*diff(y2(t), t$2) + a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)*diff(y2(t), t)^2*diff(y3(t), t) + a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)^2*diff(y3(t), t$2) + 4*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)*diff(y3(t), t)*diff(y2(t), t$2) + a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1*r_8_k1*diff(y2(t), t)^2*diff(y3(t), t)^2 + a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1*diff(y2(t), t)^2*diff(y3(t), t)*diff(y3(t), t$2) + 2*a1^2*a2^2*a3^2*r_2_k1*r_6_k1^2*r_7_k1*diff(y2(t), t)*diff(y3(t), t)^2*diff(y2(t), t$2) + 3*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)^2*y3(t)*diff(y3(t), t) + 3*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^3*r_7_k1^2*y2(t)^2*diff(y3(t), t)^2 + 2*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t) + 2*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^3*r_7_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t)^2 + 3*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y2(t)^2*y3(t)*diff(y3(t), t) + 3*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)^2*y3(t)*diff(y3(t), t$2) + 3*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y2(t)^2*diff(y3(t), t)^2 + 3*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1^2*y2(t)^2*diff(y3(t), t)*diff(y3(t), t$2) + 2*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y2(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t) + 2*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1*y2(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t$2) + 2*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t)^2 + 2*a1^2*a2^2*a3^2*r_3_k1*r_5_k1*r_6_k1^2*r_7_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t)*diff(y3(t), t$2) - a1^2*a2*a3^4*r_3_k1*r_6_k1^4*r_7_k1^2*y2(t)^2*diff(y2(t), t$2) + a1^2*a2*a3^4*r_3_k1*r_6_k1^4*r_7_k1^2*y2(t)*diff(y2(t), t)^2 - a1^2*a2*a3^4*r_3_k1*r_6_k1^4*r_7_k1*y2(t)^2*diff(y2(t), t$3) + a1^2*a2*a3^4*r_3_k1*r_6_k1^4*r_7_k1*y2(t)*diff(y2(t), t)*diff(y2(t), t$2) - a1^2*a2*a3^4*r_3_k1*r_6_k1^4*y2(t)*diff(y2(t), t)*diff(y2(t), t$3) + a1^2*a2*a3^4*r_3_k1*r_6_k1^4*y2(t)*diff(y2(t), t$2)^2 - 3*a1^2*a2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)^2*y3(t)*diff(y2(t), t) - 3*a1^2*a2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1^2*y2(t)^2*diff(y2(t), t)*diff(y3(t), t) - 4*a1^2*a2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*y3(t)*diff(y2(t), t)^2 - 4*a1^2*a2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1*y2(t)*diff(y2(t), t)^2*diff(y3(t), t) - a1^2*a2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_8_k1*y3(t)*diff(y2(t), t)^3 - a1^2*a2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*diff(y2(t), t)^3*diff(y3(t), t) - a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^3*r_8_k1*y2(t)^3*diff(y3(t), t) - 3*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^3*r_8_k1*y2(t)^2*y3(t)*diff(y2(t), t) - a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^3*y2(t)^3*diff(y3(t), t$2) - 3*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^3*y2(t)^2*diff(y2(t), t)*diff(y3(t), t) - 3*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)^2*y3(t)*diff(y2(t), t$2) - 2*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)^2*diff(y2(t), t)*diff(y3(t), t) - 4*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)*y3(t)*diff(y2(t), t)^2 - 2*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^2*y2(t)^2*diff(y2(t), t)*diff(y3(t), t$2) - 3*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^2*y2(t)^2*diff(y3(t), t)*diff(y2(t), t$2) - 4*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1^2*y2(t)*diff(y2(t), t)^2*diff(y3(t), t) - 4*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*y3(t)*diff(y2(t), t)*diff(y2(t), t$2) - a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*diff(y2(t), t)^2*diff(y3(t), t) - a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1*r_8_k1*y3(t)*diff(y2(t), t)^3 - a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1*y2(t)*diff(y2(t), t)^2*diff(y3(t), t$2) - 4*a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1*y2(t)*diff(y2(t), t)*diff(y3(t), t)*diff(y2(t), t$2) - a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_7_k1*diff(y2(t), t)^3*diff(y3(t), t) - a1^2*a2*a3^3*r_2_k1*r_6_k1^3*r_8_k1*y3(t)*diff(y2(t), t)^2*diff(y2(t), t$2) - a1^2*a2*a3^3*r_2_k1*r_6_k1^3*diff(y2(t), t)^2*diff(y3(t), t)*diff(y2(t), t$2) - a1^2*a2*a3^3*r_3_k1*r_5_k1*r_6_k1^4*r_7_k1^2*y2(t)^3*diff(y3(t), t) - a1^2*a2*a3^3*r_3_k1*r_5_k1*r_6_k1^4*r_7_k1*y2(t)^2*diff(y2(t), t)*diff(y3(t), t) - a1^2*a2*a3^3*r_3_k1*r_5_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y2(t)^3*diff(y3(t), t) - a1^2*a2*a3^3*r_3_k1*r_5_k1*r_6_k1^3*r_7_k1^2*y2(t)^3*diff(y3(t), t$2) - a1^2*a2*a3^3*r_3_k1*r_5_k1*r_6_k1^3*r_7_k1*r_8_k1*y2(t)^2*diff(y2(t), t)*diff(y3(t), t) - a1^2*a2*a3^3*r_3_k1*r_5_k1*r_6_k1^3*r_7_k1*y2(t)^2*diff(y2(t), t)*diff(y3(t), t$2) + a1^2*a3^4*r_2_k1*r_3_k1*r_6_k1^4*r_7_k1^2*y2(t)^3*diff(y2(t), t) + 2*a1^2*a3^4*r_2_k1*r_3_k1*r_6_k1^4*r_7_k1*y2(t)^2*diff(y2(t), t)^2 + a1^2*a3^4*r_2_k1*r_3_k1*r_6_k1^4*y2(t)*diff(y2(t), t)^3 + a1^2*a3^4*r_2_k1*r_6_k1^4*r_7_k1^3*y2(t)^3*diff(y2(t), t) + a1^2*a3^4*r_2_k1*r_6_k1^4*r_7_k1^2*y2(t)^3*diff(y2(t), t$2) + 2*a1^2*a3^4*r_2_k1*r_6_k1^4*r_7_k1^2*y2(t)^2*diff(y2(t), t)^2 + 2*a1^2*a3^4*r_2_k1*r_6_k1^4*r_7_k1*y2(t)^2*diff(y2(t), t)*diff(y2(t), t$2) + a1^2*a3^4*r_2_k1*r_6_k1^4*r_7_k1*y2(t)*diff(y2(t), t)^3 + a1^2*a3^4*r_2_k1*r_6_k1^4*y2(t)*diff(y2(t), t)^2*diff(y2(t), t$2) + a1*a2^4*a3^2*r_2_k2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)*y3(t)*diff(y3(t), t) - a1*a2^4*a3^2*r_2_k2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)^2*diff(y1(t), t) + a1*a2^4*a3^2*r_2_k2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y3(t), t$2) + a1*a2^4*a3^2*r_2_k2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*diff(y3(t), t)^2 - 2*a1*a2^4*a3^2*r_2_k2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y1(t), t)*diff(y3(t), t) + a1*a2^4*a3^2*r_2_k2*r_3_k1^2*r_6_k1^2*r_7_k1*y1(t)*diff(y3(t), t)*diff(y3(t), t$2) - a1*a2^4*a3^2*r_2_k2*r_3_k1^2*r_6_k1^2*r_7_k1*diff(y1(t), t)*diff(y3(t), t)^2 + a1*a2^4*a3^2*r_3_k1^3*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)*y3(t)*diff(y3(t), t) - a1*a2^4*a3^2*r_3_k1^3*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)^2*diff(y1(t), t) + a1*a2^4*a3^2*r_3_k1^3*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y3(t), t$2) + a1*a2^4*a3^2*r_3_k1^3*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*diff(y3(t), t)^2 - 2*a1*a2^4*a3^2*r_3_k1^3*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y1(t), t)*diff(y3(t), t) + a1*a2^4*a3^2*r_3_k1^3*r_6_k1^2*r_7_k1*y1(t)*diff(y3(t), t)*diff(y3(t), t$2) - a1*a2^4*a3^2*r_3_k1^3*r_6_k1^2*r_7_k1*diff(y1(t), t)*diff(y3(t), t)^2 + a1*a2^4*a3^2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)*diff(y3(t), t)^2 - a1*a2^4*a3^2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1^2*y3(t)*diff(y1(t), t)*diff(y3(t), t) + 2*a1*a2^4*a3^2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*diff(y3(t), t)*diff(y3(t), t$2) - a1*a2^4*a3^2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*y3(t)*diff(y1(t), t)*diff(y3(t), t$2) - a1*a2^4*a3^2*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*diff(y1(t), t)*diff(y3(t), t)^2 + a1*a2^4*a3^2*r_3_k1^2*r_6_k1^2*r_7_k1*y1(t)*diff(y3(t), t$2)^2 - a1*a2^4*a3^2*r_3_k1^2*r_6_k1^2*r_7_k1*diff(y1(t), t)*diff(y3(t), t)*diff(y3(t), t$2) + 2*a1*a2^4*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*r_8_k1^3*y1(t)*y3(t)^2*diff(y3(t), t) + 2*a1*a2^4*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y1(t)*y3(t)^2*diff(y3(t), t$2) + 4*a1*a2^4*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*r_8_k1^2*y1(t)*y3(t)*diff(y3(t), t)^2 + 4*a1*a2^4*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*r_8_k1*y1(t)*y3(t)*diff(y3(t), t)*diff(y3(t), t$2) + 2*a1*a2^4*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*r_8_k1*y1(t)*diff(y3(t), t)^3 + 2*a1*a2^4*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2*y1(t)*diff(y3(t), t)^2*diff(y3(t), t$2) + a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)*y3(t)^2*diff(y3(t), t) + 2*a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y3(t), t)^2 + a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1^2*r_7_k1*y1(t)*diff(y3(t), t)^3 + a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1*r_7_k1*r_8_k1^3*y1(t)*y3(t)^2*diff(y3(t), t) + a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1*r_7_k1*r_8_k1^2*y1(t)*y3(t)^2*diff(y3(t), t$2) + 2*a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1*r_7_k1*r_8_k1^2*y1(t)*y3(t)*diff(y3(t), t)^2 + 2*a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y3(t), t)*diff(y3(t), t$2) + a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1*r_7_k1*r_8_k1*y1(t)*diff(y3(t), t)^3 + a1*a2^4*a3*r_3_k1^2*r_5_k1*r_6_k1*r_7_k1*y1(t)*diff(y3(t), t)^2*diff(y3(t), t$2) - a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y2(t)*diff(y3(t), t) - a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y2(t), t) + 2*a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*y3(t)*diff(y1(t), t) - a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*r_7_k1*y1(t)*y2(t)*diff(y3(t), t$2) - a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*r_7_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t) + 2*a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*r_7_k1*y2(t)*diff(y1(t), t)*diff(y3(t), t) - a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*r_8_k1*y1(t)*y3(t)*diff(y2(t), t$2) + a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*r_8_k1*y3(t)*diff(y1(t), t)*diff(y2(t), t) - a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*y1(t)*diff(y3(t), t)*diff(y2(t), t$2) + a1*a2^3*a3^3*r_2_k2*r_3_k1^2*r_6_k1^3*diff(y1(t), t)*diff(y2(t), t)*diff(y3(t), t) - a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y2(t)*diff(y3(t), t) - a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y2(t), t) + 2*a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*y3(t)*diff(y1(t), t) - a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*r_7_k1*y1(t)*y2(t)*diff(y3(t), t$2) - a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*r_7_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t) + 2*a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*r_7_k1*y2(t)*diff(y1(t), t)*diff(y3(t), t) - a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*r_8_k1*y1(t)*y3(t)*diff(y2(t), t$2) + a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*r_8_k1*y3(t)*diff(y1(t), t)*diff(y2(t), t) - a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*y1(t)*diff(y3(t), t)*diff(y2(t), t$2) + a1*a2^3*a3^3*r_3_k1^3*r_6_k1^3*diff(y1(t), t)*diff(y2(t), t)*diff(y3(t), t) + a1*a2^3*a3^3*r_3_k1^2*r_5_k2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y2(t)*diff(y3(t), t) - a1*a2^3*a3^3*r_3_k1^2*r_5_k2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y2(t), t) + a1*a2^3*a3^3*r_3_k1^2*r_5_k2*r_6_k1^3*r_7_k1*y1(t)*y2(t)*diff(y3(t), t$2) - a1*a2^3*a3^3*r_3_k1^2*r_5_k2*r_6_k1^3*r_7_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t) + a1*a2^3*a3^3*r_3_k1^2*r_6_k1^4*r_7_k1*r_8_k1*y1(t)*y2(t)*diff(y3(t), t) - a1*a2^3*a3^3*r_3_k1^2*r_6_k1^4*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y2(t), t) + a1*a2^3*a3^3*r_3_k1^2*r_6_k1^4*r_7_k1*y1(t)*y2(t)*diff(y3(t), t$2) - a1*a2^3*a3^3*r_3_k1^2*r_6_k1^4*r_7_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t) - a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y2(t), t$2) - a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t) + a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y2(t)*diff(y1(t), t)*diff(y3(t), t) + a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y3(t)*diff(y1(t), t)*diff(y2(t), t) - a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_7_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t$2) - a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_7_k1*y1(t)*diff(y3(t), t)*diff(y2(t), t$2) + a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_7_k1*y2(t)*diff(y1(t), t)*diff(y3(t), t$2) + a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_7_k1*diff(y1(t), t)*diff(y2(t), t)*diff(y3(t), t) - a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_8_k1*y1(t)*y3(t)*diff(y2(t), t$3) + a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*r_8_k1*y3(t)*diff(y1(t), t)*diff(y2(t), t$2) - a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*y1(t)*diff(y3(t), t)*diff(y2(t), t$3) + a1*a2^3*a3^3*r_3_k1^2*r_6_k1^3*diff(y1(t), t)*diff(y3(t), t)*diff(y2(t), t$2) - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)*y3(t)^2*diff(y2(t), t) - 4*a1*a2^3*a3^2*r_2_k1*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t) - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1^2*r_6_k1^2*r_7_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t)^2 - 4*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y1(t)*y2(t)*y3(t)*diff(y3(t), t) - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1^2*y1(t)*y3(t)^2*diff(y2(t), t) - 4*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y1(t)*y2(t)*y3(t)*diff(y3(t), t$2) - 4*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y1(t)*y2(t)*diff(y3(t), t)^2 - 4*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*r_8_k1*y1(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t) - 4*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*y1(t)*y2(t)*diff(y3(t), t)*diff(y3(t), t$2) - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1^2*y1(t)*diff(y2(t), t)*diff(y3(t), t)^2 - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)*y3(t)^2*diff(y2(t), t$2) - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t) - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y2(t), t)*diff(y3(t), t$2) - 4*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y3(t), t)*diff(y2(t), t$2) - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t)^2 - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*y1(t)*diff(y2(t), t)*diff(y3(t), t)*diff(y3(t), t$2) - 2*a1*a2^3*a3^2*r_2_k1*r_3_k1*r_6_k1^2*r_7_k1*y1(t)*diff(y3(t), t)^2*diff(y2(t), t$2) - 2*a1*a2^3*a3^2*r_3_k1^2*r_5_k1*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y2(t)*y3(t)*diff(y3(t), t) - 2*a1*a2^3*a3^2*r_3_k1^2*r_5_k1*r_6_k1^3*r_7_k1*y1(t)*y2(t)*diff(y3(t), t)^2 - 2*a1*a2^3*a3^2*r_3_k1^2*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)*y2(t)*y3(t)*diff(y3(t), t) - 2*a1*a2^3*a3^2*r_3_k1^2*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*y2(t)*y3(t)*diff(y3(t), t$2) - 2*a1*a2^3*a3^2*r_3_k1^2*r_5_k1*r_6_k1^2*r_7_k1*r_8_k1*y1(t)*y2(t)*diff(y3(t), t)^2 - 2*a1*a2^3*a3^2*r_3_k1^2*r_5_k1*r_6_k1^2*r_7_k1*y1(t)*y2(t)*diff(y3(t), t)*diff(y3(t), t$2) + a1*a2^2*a3^4*r_2_k2*r_3_k1^2*r_6_k1^4*r_7_k1*y1(t)*y2(t)*diff(y2(t), t) - a1*a2^2*a3^4*r_2_k2*r_3_k1^2*r_6_k1^4*r_7_k1*y2(t)^2*diff(y1(t), t) + a1*a2^2*a3^4*r_2_k2*r_3_k1^2*r_6_k1^4*y1(t)*y2(t)*diff(y2(t), t$2) - a1*a2^2*a3^4*r_2_k2*r_3_k1^2*r_6_k1^4*y2(t)*diff(y1(t), t)*diff(y2(t), t) + a1*a2^2*a3^4*r_3_k1^3*r_6_k1^4*r_7_k1*y1(t)*y2(t)*diff(y2(t), t) - a1*a2^2*a3^4*r_3_k1^3*r_6_k1^4*r_7_k1*y2(t)^2*diff(y1(t), t) + a1*a2^2*a3^4*r_3_k1^3*r_6_k1^4*y1(t)*y2(t)*diff(y2(t), t$2) - a1*a2^2*a3^4*r_3_k1^3*r_6_k1^4*y2(t)*diff(y1(t), t)*diff(y2(t), t) + a1*a2^2*a3^4*r_3_k1^2*r_6_k1^4*r_7_k1*y1(t)*y2(t)*diff(y2(t), t$2) - a1*a2^2*a3^4*r_3_k1^2*r_6_k1^4*r_7_k1*y2(t)*diff(y1(t), t)*diff(y2(t), t) + a1*a2^2*a3^4*r_3_k1^2*r_6_k1^4*y1(t)*y2(t)*diff(y2(t), t$3) - a1*a2^2*a3^4*r_3_k1^2*r_6_k1^4*y2(t)*diff(y1(t), t)*diff(y2(t), t$2) + 4*a1*a2^2*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y2(t)*y3(t)*diff(y2(t), t) + 4*a1*a2^2*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*r_7_k1*y1(t)*y2(t)*diff(y2(t), t)*diff(y3(t), t) + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*r_8_k1*y1(t)*y3(t)*diff(y2(t), t)^2 + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*y1(t)*diff(y2(t), t)^2*diff(y3(t), t) + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y1(t)*y2(t)^2*diff(y3(t), t) + 4*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1^2*r_8_k1*y1(t)*y2(t)*y3(t)*diff(y2(t), t) + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1^2*y1(t)*y2(t)^2*diff(y3(t), t$2) + 4*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1^2*y1(t)*y2(t)*diff(y2(t), t)*diff(y3(t), t) + 4*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y2(t)*y3(t)*diff(y2(t), t$2) + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y2(t)*diff(y2(t), t)*diff(y3(t), t) + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y3(t)*diff(y2(t), t)^2 + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1*y1(t)*y2(t)*diff(y2(t), t)*diff(y3(t), t$2) + 4*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1*y1(t)*y2(t)*diff(y3(t), t)*diff(y2(t), t$2) + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_7_k1*y1(t)*diff(y2(t), t)^2*diff(y3(t), t) + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*r_8_k1*y1(t)*y3(t)*diff(y2(t), t)*diff(y2(t), t$2) + 2*a1*a2^2*a3^3*r_2_k1*r_3_k1*r_6_k1^3*y1(t)*diff(y2(t), t)*diff(y3(t), t)*diff(y2(t), t$2) + a1*a2^2*a3^3*r_3_k1^2*r_5_k1*r_6_k1^4*r_7_k1*y1(t)*y2(t)^2*diff(y3(t), t) + a1*a2^2*a3^3*r_3_k1^2*r_5_k1*r_6_k1^3*r_7_k1*r_8_k1*y1(t)*y2(t)^2*diff(y3(t), t) + a1*a2^2*a3^3*r_3_k1^2*r_5_k1*r_6_k1^3*r_7_k1*y1(t)*y2(t)^2*diff(y3(t), t$2) - 2*a1*a2*a3^4*r_2_k1*r_3_k1^2*r_6_k1^4*r_7_k1*y1(t)*y2(t)^2*diff(y2(t), t) - 2*a1*a2*a3^4*r_2_k1*r_3_k1^2*r_6_k1^4*y1(t)*y2(t)*diff(y2(t), t)^2 - 2*a1*a2*a3^4*r_2_k1*r_3_k1*r_6_k1^4*r_7_k1^2*y1(t)*y2(t)^2*diff(y2(t), t) - 2*a1*a2*a3^4*r_2_k1*r_3_k1*r_6_k1^4*r_7_k1*y1(t)*y2(t)^2*diff(y2(t), t$2) - 2*a1*a2*a3^4*r_2_k1*r_3_k1*r_6_k1^4*r_7_k1*y1(t)*y2(t)*diff(y2(t), t)^2 - 2*a1*a2*a3^4*r_2_k1*r_3_k1*r_6_k1^4*y1(t)*y2(t)*diff(y2(t), t)*diff(y2(t), t$2) + a2^4*a3^2*r_2_k1*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1^2*y1(t)^2*y3(t)*diff(y3(t), t) + a2^4*a3^2*r_2_k1*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*y1(t)^2*y3(t)*diff(y3(t), t$2) + a2^4*a3^2*r_2_k1*r_3_k1^2*r_6_k1^2*r_7_k1*r_8_k1*y1(t)^2*diff(y3(t), t)^2 + a2^4*a3^2*r_2_k1*r_3_k1^2*r_6_k1^2*r_7_k1*y1(t)^2*diff(y3(t), t)*diff(y3(t), t$2) - a2^3*a3^3*r_2_k1*r_3_k1^3*r_6_k1^3*r_8_k1*y1(t)^2*y3(t)*diff(y2(t), t) - a2^3*a3^3*r_2_k1*r_3_k1^3*r_6_k1^3*y1(t)^2*diff(y2(t), t)*diff(y3(t), t) - a2^3*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)^2*y2(t)*diff(y3(t), t) - a2^3*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*r_7_k1*r_8_k1*y1(t)^2*y3(t)*diff(y2(t), t) - a2^3*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*r_7_k1*y1(t)^2*y2(t)*diff(y3(t), t$2) - a2^3*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*r_7_k1*y1(t)^2*diff(y2(t), t)*diff(y3(t), t) - a2^3*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*r_8_k1*y1(t)^2*y3(t)*diff(y2(t), t$2) - a2^3*a3^3*r_2_k1*r_3_k1^2*r_6_k1^3*y1(t)^2*diff(y3(t), t)*diff(y2(t), t$2) + a2^2*a3^4*r_2_k1*r_3_k1^3*r_6_k1^4*y1(t)^2*y2(t)*diff(y2(t), t) + a2^2*a3^4*r_2_k1*r_3_k1^2*r_6_k1^4*r_7_k1*y1(t)^2*y2(t)*diff(y2(t), t) + a2^2*a3^4*r_2_k1*r_3_k1^2*r_6_k1^4*y1(t)^2*y2(t)*diff(y2(t), t$2)]:
#IOeqs := [expand(eqs[1]/(a2*a3^2*r_6_k1^2)), expand(eqs[2]/(a2*a3*r_3_k1*r_6_k1)), expand(eqs[3]/(2*a1*a2^4*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2))]:
IOeqs := eqs:

#Replacing derivatives with regular variables for further computation

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

IOeqs_norm := [expand(eqs[1]/(a2*a3^2*r_6_k1^2)), expand(eqs[2]/(a2*a3*r_3_k1*r_6_k1)), expand(eqs[3]/(2*a1*a2^4*a3*r_2_k1*r_3_k1*r_6_k1*r_7_k1^2))]:

IO_coeffs := [seq(i, i in coeffs(expand(IOeqs_norm[1]), map(lhs, alg_subs))), seq(i, i in coeffs(expand(IOeqs_norm[2]), map(lhs, alg_subs))), seq(i, i in coeffs(expand(IOeqs_norm[3]), map(lhs, alg_subs)))]:
IO_coeffs := FilterGenerators(FieldToIdeal(IO_coeffs)):
IO_coeffs_nonumbers := []:

for i in IO_coeffs do 
  if not type(i, integer) then IO_coeffs_nonumbers := [op(IO_coeffs_nonumbers), i] fi; 
od:

for i in IO_coeffs_nonumbers do 
  printf("%a\n", i) 
od:
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


printf("\n======================================\n");
printf("Computation for choosing alpha-tildes:\n");

ov := Flatten(LieDer):
LieDerNo_t := subs([seq(xv(t) = xv, xv in [op(states), op(inputs)])], ov):
Jac := simplify(Jacobian(LieDerNo_t, states)):
size_minor := add(degree(outputs_maxorders[num_out][1]), num_out = 1..numelems(outputs)):
Q := Determinant(Jac[1..size_minor, 2..size_minor+1]):
   
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
tilde_eqs := [seq(subs(params_sub, i) = i, i in IO_coeffs_nonumbers)]:
tilde_system := [op(tilde_eqs), subs(params_sub, Q0) <> 0]:

printf("\nThe corresponding system to choose the tilde-parameters is:\n\n");

for s in tilde_system do printf("%a\n", s) od;

printf("\nWe pick the following solution:\n\n");

chosen_tildes := solve([op(tilde_system), 
                        X1_t_tilde = 1, r_9_k1_tilde=1, r_5_k1_tilde = 1, r_1_k2_tilde = -1], params_tilde)[1]:

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
for i from 1 to numelems(LieDer) do
  for j from 1 to numelems(LieDer[i]) do 
   eq_conv[i, j] := LieDer_no_t[i][j] = LieDerW[i][j]
  od; 
od;
conv_sols := simplify(solve([seq(seq(eq_conv[i, j], j=1..numelems(LieDer[i])), i = 1..numelems(LieDer))], wvars)):
printf("\n=========================================\n");
printf("The corresponding change of variables is:\n\n");

for s in conv_sols[1] do printf("%a\n", s): od;



