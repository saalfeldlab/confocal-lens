# confocal-lens
This repository contains some Java classes and scripts to calculate lens-distortion models including chromatic abberation from confocal stacks and how to apply them.

If you are using this tool in your work, please do not forget to cite [our paper](https://www.biorxiv.org/content/10.1101/376384v1) where we first introduce this method:

JA Bogovic, H Otsuna, L Heinrich, M Ito, J Jeter, GW Meissner, A Nern, J Colonell, O Malkesman, K Ito, S Saalfeld, "An unbiased template of the Drosophila brain and ventral nerve cord", bioRxiv 376384, DOI: [10.1101/376384](https://doi.org/10.1101/376384).

# Multi-color lens distortion and shift correction
Please follow the instructions in this [screencast](https://www.youtube.com/watch?v=lPt-WQuniUs).  The approach is equivalent to what we have done for multi-camera calibration across electron microscopes.

# Apply distortion models

This is an example for the example models included in this repository.  Please modify these models and the scripts according to your specific needs.

Download the [json files](https://github.com/saalfeldlab/confocal-lens/tree/master/scripts) (multi-channel lens-models for all scopes) and drop the file [apply-lens.bsh](https://github.com/saalfeldlab/confocal-lens/blob/master/scripts/apply-lens.bsh) into the plugin directory of your Fiji installation.  After a restart of Fiji, you should find a menu-entry for this script in the plugins menu.  Click it.

If your Fiji is not too old, then it should open a dialog that asks you
for two input image paths (expected is

1. LSM file with two channels (488 and 594nm)
2. LSM file with three channels (488, 561, and 647nm)

), the respective JSON file for the scope that was used, and an output directory to save the resulting multi-channel file.  The output will have all five channels in a single file.

The script currently assumes that the lens-models in the JSON file and the channels from input images are in the same order (in this case 488, 594, 488, 561, 647nm).  You can do subsets or other combinations but you would have to edit the JSON file accordingly.  The JSON files are simple to read, check

https://github.com/saalfeldlab/confocal-lens/blob/master/scripts/scope1.json

and they contain name tags for each model.  That should make editing them easier.  The script ignores the name tags and goes by the order only.  The script should be macro-recordable.
