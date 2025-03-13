Test dataset for training of models mapping charge to electronegativity.

To generate the data by DFT simulation:

    generate_reference_data.sh
    run_cp2k.sh

Then run the generated list of CP2k jobs.

For analysis of the charges and energies, separated into a slow parsing step and a fast analysis step:

    extract_charge_energy.sh
    analyze_cp2k_charge.sh

And extra analysis of the cube files, which is slower and not needed for analysis of the charges:

    analyze_cp2k_cube.sh

The environment variables that need to be set are

    ELECTRONEGATIVITYSTRUCTUREPATH
    CP2KPSEUDOPOTENTIALPATH
    CP2KCOMMAND
