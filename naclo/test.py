import pandas as pd


a = pd.read_excel('Retrieve_smiles.xlsx')
b = pd.read_excel('6_23_CPI_lib_no_barcode_patricia.sdf.xlsx')
print(a)

substances = list(a.SUBSTANCE_ID)


substances = ['CPI00' + i[3:] for i in substances]
print(substances)
print(len(substances))

sery = b['SUBSTANCE_ID']
b2 = b.mask(~sery.isin(substances)).dropna()
print(len(b2))

b2.to_excel('retrieved_smiles_out.xlsx')



# other_isin = ref_data[globe.compare_column_name].isin(oth_data[globe.compare_column_name])
#         isin_series = pd.Series(other_isin)