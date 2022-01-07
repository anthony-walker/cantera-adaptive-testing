# zero threshold
export PRECON_OPTS="$CURR_MODEL -L -v -M -P -S GMRES -T 0 $ADD_ARGS"
echo $PRECON_OPTS
# varying thresholds
for th in {1..18}
do
    echo $PRECON_OPTS
    export PRECON_OPTS="$CURR_MODEL -L -v -M -P -S GMRES -T 1e-$th $ADD_ARGS"
    sleep 0.25
done