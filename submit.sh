# seed=31

# length=50

iterations=10000000
# dt=0.01

kB=1
# T=1

# D=0.017

LJ_e=0.1
LJ_sigma=0.5
LJ_rc=4.0

bond_k=100
bond_r0=1.0

steps_to_skip=1000

bond_break_message=0
stop_on_breakage=0

input=input.inp
output=out.run

###END OF FIXED PARAMETERS###


for length in 50
do #Enter length

mkdir -p len_${length}
cd len_${length}

for seed in {1..20}
do #Enter seed

mkdir -p seed_${seed}
cd seed_${seed}

for dt in 0.01
do #Enter dt

mkdir -p dt_${dt}
cd dt_${dt}

for T in 1.0
do #Enter T

mkdir -p T_${T}
cd T_${T}

for D in 0.017
do #Enter D

mkdir -p D_${D}
cd D_${D}

rm -f plysim
cp ../../../../../polysim .

rm -f ${input}
cp ../../../../../template_input.inp ${input}

sed -i -e "s/VAL_SEED/${seed}/g" ${input}
sed -i -e "s/VAL_LENGTH/${length}/g" ${input}
sed -i -e "s/VAL_ITER/${iterations}/g" ${input}
sed -i -e "s/VAL_DT/${dt}/g" ${input}
sed -i -e "s/VAL_KB/${kB}/g" ${input}
sed -i -e "s/VAL_T/${T}/g" ${input}
sed -i -e "s/VAL_D/${D}/g" ${input}
sed -i -e "s/VAL_LJE/${LJ_e}/g" ${input}
sed -i -e "s/VAL_LJSIGMA/${LJ_sigma}/g" ${input}
sed -i -e "s/VAL_RC/${LJ_rc}/g" ${input}
sed -i -e "s/VAL_BONDK/${bond_k}/g" ${input}
sed -i -e "s/VAL_BONDRZERO/${bond_r0}/g" ${input}
sed -i -e "s/VAL_STEPSTOSKIP/${steps_to_skip}/g" ${input}
sed -i -e "s/VAL_BREAKMSG/${bond_break_message}/g" ${input}
sed -i -e "s/VAL_STOPONBRK/${stop_on_breakage}/g" ${input}

echo Running L=${length} seed=${seed} dt=${dt} T=${T} D=${D}

# time ./polysim ${input} > ${output}

cd ../
done #Return to T

cd ../
done #Return to dt

cd ../
done #Return to seed

cd ../
done #Return to length

cd ../
done #Return to top level
