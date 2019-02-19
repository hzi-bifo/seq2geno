FILES=$(grep -rn '../../smk' -e 'source'| grep -v 'old\/'| grep ' activate' | \
  awk -F':' '{print $1}')
for f in $FILES; do
  new_f=$(echo $f| awk -F'\/' '{print $NF}')
  echo $new_f
  cat $f
#  cat $f| sed '/source activate /s/source activate /source activate $SEQ2GENO_HOME\/env\//' > ./$new_f
#  mv ./$new_f $f
done
