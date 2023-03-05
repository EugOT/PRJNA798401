npr <- c(
    "Oxtr", "Avpr1a", "Avpr1b", "Avpr2",
    "Tacr1", "Mc1r", "Mc3r", "Mc4r", "Oprl1", "Tacr3",
    "Gpr83", "Galr1", "Npy1r", "Npy2r", "Npy4r", "Npy4r2", "Npy5r",
    "Sstr1", "Sstr2", "Sstr3", "Mchr1", "Oprd1", "Oprk1", "Oprm1",
    "Trhr", "Hcrtr1", "Hcrtr2", "Qrfpr", "Npffr1", "Npffr2",
    "Prlhr", "Ghr", "Ghrhr", "Grpr", "Vipr1", "Vipr2", "Prokr2",
    "Nmur2", "Nmur1", "Nmbr", "Kiss1r", "Crhr1", "Crhr2",
    "Cntfr", "Cckar", "Cckbr", "Galr3"
)

np <- c(
    "Adcyap1", "Oxt", "Avp", "Tac1", "Pomc", "Pnoc", "Tac2",
    "Nts", "Gal", "Agrp", "Npy", "Sst", "Cartpt", "Pmch",
    "Reln", "Rxfp1", "Penk", "Pdyn", "Trh", "Hcrt", "Qrfp",
    "Npw", "Npvf", "Ghrh", "Grp", "Vip", "Nms", "Nmu", "Nmb",
    "Kiss1", "Crh", "Bdnf", "Cntf", "Cck"
)

irs_genes <- c(
    "Alk", "Insr", "Ltk", "Igf1r", "Irs1",
    "Ptn", "Mdk", "Fam150a", "Fam150b",
    "Mc4r", "Lepr", "Sim1", "Lmo4",
    "Slc2a1", "Slc2a3"
)

neurotrans <- c(
    "Slc17a6", "Slc17a7", "Slc17a8", "Slc1a1", "Slc1a2", "Slc1a6",
    "Gad1", "Slc32a1", "Slc6a1"
)
glut <- c("Slc17a6", "Slc17a7", "Slc17a8", "Slc1a1", "Slc1a2", "Slc1a6")
gaba <- c("Gad1", "Gad2", "Slc32a1", "Slc6a1")
dopam <- c("Th", "Slc6a3", "Slc18a2", "Ddc", "Slc18a3")
ach <- c("Chat", "Slc18a3", "Ache", "Slc5a7")
mcr_genes <- c("Mc1r", "Mc2r", "Mc3r", "Mc4r")
pvn_genes <- c("Trh", "Crh", "Mbnl3", "Pgf", "Irs4", "Gpr101", "Nr3c2", "Agtr1")
dmh_trh_g <- c("Onecut3", "Cartpt")

genes.embed <- c(
    "Npy1r", "Npy2r", "Esr1", "Prlr", "Otp", "Adcyap1r1", "Htra1",
    "Mc4r", "Mc3r", "Mc2r", "Mc1r", "Lxn", "Mfn2", "Gfap", "Nfia", "Fgfr3",
    "Cst3", "Slc38a1", "Ntsr2", "Apoe", "Slit2", "Ndrg2", "Plcb1", "Egfr",
    "Aldh1l1", "Aldh1a1", "Tafa1", "Sgcd", "Slc1a3", "S100b", "S100a6",
    "Ntrk2", "Glul", "Aldoc", "Gja1", "Agt", "Slc6a11", "Aqp4", "Sox9"
)

genes.nature <- c(
    "Ndrg2", "Aqp4", "Gja1",
    "Slc17a9", "Slc17a7", "Slc17a6",
    "Slc16a1", "Slc16a7", "Slc16a3",
    "Panx1", "P2rx7", "Srr",
    "Dao", "Vamp2", "Gab1",
    "Slc18a1", "Slc18a2", "Slc18a3",
    "Slc17a5", "Osmr", "S100a6",
    "Ogn", "Itih5", "Rdh10",
    "Fst", "1500015O10Rik", "Rnf13",
    "A2m", "Rsph9", "Galnt16",
    "Rad23b", "Tgfbr2", "Ppp1r15a",
    "Mlc1", "Slc6a11", "Slc1a3", "S100b",
    "S100b", "Fabp7", "Gab1", "Igfbp5",
    "Hopx", "Igsf1", "Tgfb2", "2810459M11Rik",
    "Rnf13", "Itih5", "Slc1a2", "Cd59a",
    "Vim", "Slc7a10", "Fos"
)

genes.jj <- c(
    "Acsbg1", "Acsl3", "Actb", "Agt", "Aldh1l1", "Aldoc", "Apoe",
    "Aqp4", "Atp1a2", "Bmpr1b", "Cbs", "Ckb", "Clu", "Cnx43",
    "Cpe", "Cst3", "Cth", "Cyp4f14", "Dbi", "Dbx2", "Dctd",
    "Dio2", "Ednrb", "Fgfr3", "Gabbr1", "Gabbr2", "Gjb6", "Gli3",
    "Gm266", "Gpr37l1", "Grhl1", "Grm3", "Gs", "Gucy2c", "Hapln1",
    "Hes5", "Heyl", "Hgf", "Id2", "Itih3", "Lars2", "Lcat",
    "Ldhb", "Malat1", "Mc3r", "Mfge8", "Mlc1", "Mmd2", "Mt1",
    "Mt2", "Ndrg2", "Nfia", "Npas3", "Nrxn1", "Ntrk2", "Ntsr2",
    "Olx1", "Otx2", "Pax6", "Pbxip1", "Phkg1", "Pla2g3", "Pla2g7",
    "Plcd4", "Plce1", "Plp1", "Ppap2b", "Ppia", "Prodh", "Ptprz1",
    "Rfx4", "Rpl41", "S1pr1", "Scd2", "Serpinb1c", "Slc19a3", "Slc1a2",
    "Slc1a3", "Slc39a12", "Slc4a4", "Slc6a11", "Son", "Sox2", "Sox9",
    "Sparc", "Sparcl1", "Tlr3"
)

genes.anatomy.jj <- c(
    "Gfap", "Meg3", "Snhg11", "Ttc3", "Cst3", "Sparcl1",
    "Ndrg2", "Cnx43", "Agt", "S100b", "Aldh1l1", "Lxn",
    "Aqp4", "Slc1a2", "Fabp7", "Glul", "Olig2"
)

# public resources:
housekeeping_mouse <-
    read_lines(file = here(data_dir, "housekeeping_mouse.tsv"))
astroenriched_mouse <-
    read_lines(file = here(data_dir, "astrocyte_enriched_genes_shared_by_mice.tsv"))

astroprogenitor_humans <-
    read_lines(file = here(data_dir, "top_astrocytes_progenitor_cells_genes_by_humans.tsv")) |>
    str_to_sentence()
astromature_humans <-
    read_lines(file = here(data_dir, "top_mature_astrocyte_genes_by_humans.tsv")) |>
    str_to_sentence()

# aggregate:
genes.manual <- unique(c(
    genes.embed, genes.nature,
    genes.jj, genes.anatomy.jj,
    astroenriched_mouse,
    astromature_humans,
    astroprogenitor_humans
))


gene_int <-
    c(
        npr, np, irs_genes, neurotrans, mcr_genes,
        genes.manual, pvn_genes, dmh_trh_g
    ) %>% unique()
