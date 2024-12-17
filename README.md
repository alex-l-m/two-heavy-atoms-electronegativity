Test dataset for training of models mapping charge to electronegativity.

To generate the data by DFT simulation:

    generate_reference_data.sh
    run_cp2k.sh

Then run the generated list of CP2k jobs.

For analysis of the data:

    analyze_cp2k_cube.sh
    parse_cp2k_logs.sh
    analyze_cp2k_charge.sh

The environment variables that need to be set are

    ELECTRONEGATIVITYSTRUCTUREPATH
    CP2KPSEUDOPOTENTIALPATH
    CP2KCOMMAND
