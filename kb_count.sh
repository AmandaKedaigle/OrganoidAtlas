# Kallisto bustools count 
# specifically for RNA velocity in the downstream

ref_dir=$1
outname=$2
read_dir=$3

echo ${outname}
echo ${read_dir}

reads=$(ls ${read_dir}/*R*)

echo ${reads}


kb count \
--h5ad \
-i ${ref_dir}/index.idx \
-g ${ref_dir}/t2g.txt \
-x 10xv3 \
-w 3M-february-2018.txt \
-o ../kb_counts/${outname} \
-c1 ${ref_dir}/cdna_t2c.txt \
-c2 ${ref_dir}/intron_t2c.txt \
--workflow lamanno \
${reads}


echo DONE

