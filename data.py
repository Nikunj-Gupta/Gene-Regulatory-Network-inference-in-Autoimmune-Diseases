import GEOparse, pandas as pd, os 



# gse = GEOparse.get_GEO(geo="GSE1563", destdir="./")
gse = GEOparse.get_GEO(filepath="./datasets/GSE32591_family.soft.gz")


# print()
# print("GSM example:")
# for gsm_name, gsm in gse.gsms.items():
#     print("Name: ", gsm_name)
#     print("Metadata:",)
#     for key, value in gsm.metadata.items():
#         print(" - %s : %s" % (key, ", ".join(value)))
#     print ("Table data:",)
#     print (gsm.table.head())
#     break 

print()
print("GPL example:")
for gpl_name, gpl in gse.gpls.items():
    print("Name: ", gpl_name)
    print("Metadata:",)
    for key, value in gpl.metadata.items():
        print(" - %s : %s" % (key, ", ".join(value)))
    print("Table data:",)
    print(gpl.table) 
    print(gpl.table.columns)
    break


# df_tag = False 
# for gsm_name, gsm in gse.gsms.items():
#     # df = gsm.table 
#     # print(df[df['ID_REF'].str.contains('9145', regex=False)]) 
#     if not df_tag: 
#         df = gsm.table 
#         df_tag = True 
#     else: 
#         df = pd.merge(df, gsm.table, on='ID_REF') 


#     # break 
# print(df) 


# import cdflib 
# cdf_file = cdflib.CDF('./datasets/GSE32591_RAW/GPL14663_Hs133A_Hs_ENTREZG.cdf')
# print(cdf_file.cdf_info())