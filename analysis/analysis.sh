object_name=$( (rev | cut -d / -f 2 | rev) <<< "$PWD")

python roc.py "$object_name" highpass
python roc.py "$object_name" smooth
python roc.py "$object_name" numbasis
python maxheatmap.py "$object_name" annuli subsections
python maxheatmap.py "$object_name" annuli movement
python maxheatmap.py "$object_name" subsections movement
python meanheatmap.py "$object_name" annuli subsections
python meanheatmap.py "$object_name" annuli movement
python meanheatmap.py "$object_name" subsections movement
