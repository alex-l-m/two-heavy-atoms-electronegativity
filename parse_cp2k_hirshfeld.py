'''Read Hirshfeld charges from CP2K output'''
import re
from os.path import basename, splitext, exists
import pandas as pd

charge_section_regex = re.compile(r''' *!-+! *
 *Hirshfeld +Charges *
 *
 *#Atom +Element +Kind +Ref +Charge +Population +Net +charge *
 *1 +([A-Z][a-z]?) +1 +[-+.0-9]+ +[-+.0-9]+ +([-+.0-9]+) *
 *2 +([A-Z][a-z]?) +2 +[-+.0-9]+ +[-+.0-9]+ +([-+.0-9]+) *
 *
 *Total +Charge +[-+.0-9]* *
 *!-+! *''')

infiles = ['cp2k_logs/zincblende_Ga_As_F0_Vfield.out']

intable = pd.read_csv('simulations.csv')

# Rows of the table as a list of dictionaries
outrows = []
for row in intable.itertuples():
    simulation_id = row.simulation_id
    logfile_path = row.log_file_path
    # Skip if the file doesn't exist
    if not exists(logfile_path):
        continue
    simulation_id, ext = splitext(basename(logfile_path))
    logfile_text = open(logfile_path).read()
    match = re.search(charge_section_regex, logfile_text)
    donor_element, donor_charge_str, \
        acceptor_element, acceptor_charge_str = match.groups()
    outrows.append({'simulation_id': simulation_id, 'symbol': donor_element,
                 'donor_or_acceptor': 'donor',
                 'charge': donor_charge_str})
    outrows.append({'simulation_id': simulation_id, 'symbol': acceptor_element,
                 'donor_or_acceptor': 'acceptor',
                 'charge': acceptor_charge_str})

outtbl = pd.DataFrame(outrows)
outtbl.to_csv('cp2k_hirshfeld_charges.csv.gz', index=False)
