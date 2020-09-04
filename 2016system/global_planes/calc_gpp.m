%, , cphiNH, cthetaNH, cphiCCa, cthetaCCa, cphiCH, cthetaCH, cphiCN, cthetaCN, cphiNCA, cthetaNCA
[H N CA C O CB CG SD CE CD NE2 OE1 CD1 CD2 CE1 CE2 CZ OH NZ CG1 CG2 ND2 OD1 OG1 OE2 OD2 CE3 NE1 CZ2 CZ3 CH2 OXT ] = atom_positions()
[cphiCcsax, cthetaCcsax, cphiCcsay, cthetaCcsay, cphiNcsax, cthetaNcsax, cphiNcsay, cthetaNcsay, cphiNH, cthetaNH, cphiCCa, cthetaCCa, cphiCH, cthetaCH, cphiCN, cthetaCN, cphiNCA, cthetaNCA] = coords_global_peptide_planes(H, N, CA, C, O)

output_file(cphiCcsax, cthetaCcsax, "GPP_CCSAx")
output_file(cphiCcsay, cthetaCcsay, "GPP_CCSAy")
output_file(cphiNcsax, cthetaNcsax, "GPP_NCSAx")
output_file(cphiNcsay, cthetaNcsay, "GPP_NCSAy")
output_file(cphiNH, cthetaNH, "GPP_NH")
output_file(cphiCCa, cthetaCCa, "GPP_CCA")
output_file(cphiCH, cthetaCH, "GPP_CH")
output_file(cphiCN, cthetaCN, "GPP_CN")
output_file(cphiNCA, cthetaNCA, "GPP_NCA")


