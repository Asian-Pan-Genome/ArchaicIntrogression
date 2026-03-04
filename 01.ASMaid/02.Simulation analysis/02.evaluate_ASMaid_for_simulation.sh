#!/usr/bin/sh

out_dir="Sim_Result"
sample_dir=$out_dir/sample

function merge_all_sample_bed {
    suffix=$1
    merge_bed="$out_dir/note.merge.$suffix"
    echo "--------------------------------------"
    echo "Process $suffix ..."

    if [ -f $merge_bed ]; then
        echo "File $merge_bed already exists. Remove..."
        rm $merge_bed
    fi

    for bed in $(ls $sample_dir/*.$suffix)
    do
        sample_name=$(basename $bed | cut -d'.' -f1)
        # Substitute the first column ("chro" info) with sample name
        cat $bed | grep -v '^#' | awk -v a=$sample_name -v OFS="\t" '{ $1 = a; print }' >> $merge_bed
    done
}

merge_all_sample_bed "intro"
merge_all_sample_bed "AFR-10.tracts"
merge_all_sample_bed "AFR-10.gt005.tracts"
merge_all_sample_bed "AFR-10.gt010.tracts"
merge_all_sample_bed "AFR-10.gt015.tracts"
merge_all_sample_bed "AFR-5.tracts"
merge_all_sample_bed "AFR-15.tracts"
merge_all_sample_bed "AFR-20.tracts"


# Get Bed for different Probability
cat $out_dir/note.merge.AFR-10.tracts | awk '{if($4>=0.7) print $0}' > $out_dir/note.merge.AFR-10.prob07.tracts
cat $out_dir/note.merge.AFR-10.tracts | awk '{if($4>=0.8) print $0}' > $out_dir/note.merge.AFR-10.prob08.tracts
cat $out_dir/note.merge.AFR-10.tracts | awk '{if($4>=0.9) print $0}' > $out_dir/note.merge.AFR-10.prob09.tracts


eval_name=$out_dir/intro.eval
if [ -f $eval_name ]; then
    echo "File $eval_name already exists. Remove..."
    rm $eval_name
fi

function get_eval_ratio {
    local truth_bed=$1
    local infer_bed=$2
    
    local truth_num=$(cat $truth_bed | wc -l)
    local infer_num=$(cat $infer_bed | wc -l)
    local truth_len=$(cat $truth_bed | awk '{print $3-$2+1}' | awk '{sum+=$1} END {print sum}')
    local infer_len=$(cat $infer_bed | awk '{print $3-$2+1}' | awk '{sum+=$1} END {print sum}')

    local overlap_len=$(bedtools intersect -wao -a $truth_bed -b $infer_bed | cut -f 14 | awk '{sum+=$1} END {print sum}')
    local overlap_truth_num=$(bedtools intersect -wa -a $truth_bed -b $infer_bed | uniq | wc -l)
    local overlap_infer_num=$(bedtools intersect -wa -a $infer_bed -b $truth_bed | uniq | wc -l)
    local false_positive_bp=$(bedtools intersect -v -a $infer_bed -b $truth_bed | awk '{sum+=$3-$2+1} END {print sum}')

    local precision_len=$(echo "scale=4; $overlap_len/$infer_len" | bc)
    local precision_count=$(echo "scale=4; $overlap_infer_num/$infer_num" | bc)
    local recall_len=$(echo "scale=4; $overlap_len/$truth_len" | bc)
    local recall_count=$(echo "scale=4; $overlap_truth_num/$truth_num" | bc)
    local f1_len=$(echo "scale=4; 2*$precision_len*$recall_len/($precision_len+$recall_len)" | bc)
    local f1_count=$(echo "scale=4; 2*$precision_count*$recall_count/($precision_count+$recall_count)" | bc)
    local FPR=$(echo "scale=4;  $false_positive_bp/(100000000*100-$truth_len)" | bc)
    
    # echo "percision_len: $precision_len"
    # echo "recall_len: $recall_len"
    # echo "f1_len: $f1_len"
    # echo "percision_count: $precision_count"
    # echo "recall_count: $recall_count"
    # echo "f1_count: $f1_count"
    # echo "FPR: $FPR"
    echo -e -n "$precision_len\t$recall_len\t$f1_len\t$FPR\t$precision_count\t$recall_count\t$f1_count\t"
}

function multi_len_filter_evalution {
    local truth_bed=$1
    local infer_bed=$2
    infer_bed_label=$(basename $infer_bed | sed 's/note.merge.//g' | sed 's/.tracts//g')
    echo ">>> Process $infer_bed_label ..."
    #echo "### "$infer_bed_name >> $eval_name

    # Length filter threshold
    local lengths=(0 10000 20000 30000 40000 50000 60000)

    for len in "${lengths[@]}";do
        #cat $truth_bed | awk -v len=$len '{if($3-$2+1>=len) print $0}' > $truth_bed.$len
        cat $infer_bed | awk -v len=$len '{if($3-$2+1>=len) print $0}' > $infer_bed.$len
        echo -e -n "$len\t" >> $eval_name
        get_eval_ratio "$truth_bed" "$infer_bed.$len"  >> $eval_name
        echo $infer_bed_label >> $eval_name
        rm $infer_bed.$len
    done
}

multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-10.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-10.gt005.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-10.gt010.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-10.gt015.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-10.prob07.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-10.prob08.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-10.prob09.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-5.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-15.tracts"
multi_len_filter_evalution "$out_dir/note.merge.intro" "$out_dir/note.merge.AFR-20.tracts"