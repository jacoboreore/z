clear
echo "JACOB's fight club compiler"
echo "Version control:"

echo "2025"
echo "09/02 parametrizable settings"
echo "09/02 small z separated calculation, issues calculating small t (major change)"
echo "09/05 thread z calculations RS issues, code cleansing"
echo "09/10 Valley scanner first beta, v0.1.0, stable (milestone)"
echo "09/10 Smart step to increase walking stride"
echo "09/11 RS + Gabcke correction, single thread, R correction issues (major change)"
echo "09/11 Newton RS + Gabcke correction, single thread and threaded Newton refinement (major change)"
echo "09/13 Smart walk, better approximation to the valley"
echo "09/15 Adaptative valley walking, recursive valley approach to refine"
echo "09/16 Odlyzko–Schönhage implementation, valley walk adapatation (major change)"
echo "09/17 Plaform precision adaptation for Odlyzko–Schönhage implementation, to maintain speed and 'long double' calculations"
echo "09/17 Crashing handling, try catch handling"
echo "09/17 Stable Odlyzko–Schönhage valley scanner algorithm merge, valley scanner using pure OS since RS-Gabcke is taxing.  Typos in comments"
echo "09/18 Docker version, native compilation, CPUs as param (major adding), run.sh"
echo "09/21 Birthday :D"
echo "09/22 Q&A, bad spacing calcs, refining valley walking"
echo "09/25 Q&A, walking through valleys, stable, some zeros are not beign calculated"
echo "09/26 Q&A, fixes"
echo "09/29 Constant "step" for valley walking, this ensures to find all zeros"
echo "09/29 Constant step formula sqrt(N) / 10"

echo "-------------------------------------------"

# native compilation
make clean && make