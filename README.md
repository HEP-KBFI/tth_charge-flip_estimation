# tth_charge-flip_estimation

## Installation
See https://github.com/HEP-KBFI/tth-htt/wiki/Setup#for-auxiliary-measurements

## Running

1. Create datacards for data and pseudodata:
```bash
create_pseudodata_datacard.py -i prepareDatacards_charge_flip_mass_ll.root \
  -o prepareDatacards_data_charge_flip_mass_ll.root \
  -O prepareDatacards_pseudodata_charge_flip_mass_ll.root \
  -t ele -s data_obs DY DY_fake Singletop Diboson TTbar
```
2. Create text files from both datacards:
```bash
ChargeFlipDC -i prepareDatacards_data_charge_flip_mass_ll.root \
  -b DY_fake Singletop Diboson TTbar \
  -o fit_results_data -y 2016 -e -f

ChargeFlipDC -i prepareDatacards_pseudodata_charge_flip_mass_ll.root \
  -b DY_fake Singletop Diboson TTbar \
  -o fit_results_pseudodata -y 2016 -e -f
```
3. Fit the 21 electron pair categories and make plots:
```bash
make_fits.py -i fit_results_data       -o fit_results_data/results_cat.txt
make_fits.py -i fit_results_pseudodata -o fit_results_pseudodata/results_cat.txt
```
4. Create prefit and postfit plots:
```bash
make_fit_plots.py -i fit_results_data
make_fit_plots.py -i fit_results_pseudodata
```
5. To obtain flip rates and create various pull plots:
```bash
plot_pulls.py -f fit_results_data -F fit_results_pseudodata -i hadd_stage2_Tight.root -o $PWD -l
```

Use `fit_result_data_exclusions.root` in the analysis to assign charge misidentification probabilities to the data events.
