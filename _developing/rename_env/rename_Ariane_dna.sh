FILES=$(grep -rn '../../smk' -e 'source'| grep -v 'old\/'| grep ' activate' | \
  grep Ariane| awk -F':' '{print $1}')
for f in $FILES; do
  new_f=$(echo $f| awk -F'\/' '{print $NF}')
  echo $new_f
  cat $f| sed '/source activate Ariane_dna/s/Ariane_dna/old_mapping/' > ./$new_f
  mv ./$new_f $f
done
