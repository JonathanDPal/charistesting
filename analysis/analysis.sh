object_name=$1
num_injections=$2

python roc.py "$object_name" highpass "$num_injections"
python roc.py "$object_name" smooth "$num_injections"
python roc.py "$object_name" numbasis "$num_injections"
python maxheatmap.py "$object_name" annuli subsections
python maxheatmap.py "$object_name" annuli movement
python maxheatmap.py "$object_name" subsections movement
python meanheatmap.py "$object_name" annuli subsections "$num_injections"
python meanheatmap.py "$object_name" annuli movement "$num_injections"
python meanheatmap.py "$object_name" subsections movement "$num_injections"
