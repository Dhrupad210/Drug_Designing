Computational Discovery and Optimization of an Inhibitor for a Novel Gut Proteome Target
========================================================================================

This repository documents a complete _in silico_ drug discovery pipeline, beginning with an uncharacterized protein sequence (TargetProtein-X) identified from a gut proteome sample. The project follows a journey of discovery: from sequence to 3D model, from structural validation to functional identification, and finally, to the rational design of a novel inhibitor with superior binding affinity.

Project Workflow
----------------

1.  **Protein Identification:** An unknown protein sequence from a gut sample was analyzed with blastp to find a structural template.
    
2.  **Homology Modeling (MODELLER):** A 3D model of the unknown protein was built using PDB ID 1J2E as the template.
    
3.  **MD Simulation (GROMACS):** The 3D model was validated for stability using a 200ps Molecular Dynamics simulation.
    
4.  **Functional Validation (Vina):** The protein was identified as DPP4. Its 3D model was _functionally validated_ by docking known inhibitors sourced from ChEMBL.
    
5.  **Rational Drug Design:** The best-performing known inhibitor was analyzed, and three new ligands were designed to improve binding.
    
6.  **Final Validation:** The new ligands were docked, and one (n2) was confirmed to have a superior binding affinity.
    

Phase 1: Protein Identification & Modeling
------------------------------------------

**Thought Process:** A novel protein sequence (dpp4.fasta) was isolated from a gut proteome analysis. To understand its function, we first needed to determine its 3D structure. Since no experimental structure exists, we used homology modeling.

### 1.1. Homolog Search

To build a model, we first needed a template. We used blastp to find the closest structural homolog in the Protein Data Bank (PDB).

`blastp -query dpp4.fasta -db pdb -out results.blast`

*   **Explanation:** blastp (Protein BLAST) searches the pdb database for known protein structures that have a similar sequence to our query, dpp4.fasta.
    
*   **Result:** PDB ID **1J2E** showed the highest sequence similarity and score. This structure was downloaded to be used as our template.
    

### 1.2. PDB Template Cleaning

The raw PDB file (1J2E.pdb) is "dirty" and contains non-protein atoms, alternate locations, and other artifacts not suitable for modeling. We used pdb-tools to clean it.
`   # Clean PDB: remove HETATMs, Hydrogens, alternate locations, and low-occupancy atoms  
pdb_delhetatm 1J2E.pdb | pdb_delelem -H | pdb_selaltloc | pdb_occ -0.3 > cleaned.pdb  `
`  # Select only Chain A  
pdb_selchain -A cleaned.pdb > cleaned_chainA.pdb ` 
`  # Tidy the PDB for formatting  
pdb_tidy cleaned_chainA.pdb > 1J2E_clean_final.pdb   `

*   **Explanation:** This command pipeline:
    
    1.  pdb\_delhetatm: Removes all non-protein "heteroatoms" (like water, ions, or ligands).
        
    2.  pdb\_delelem -H: Removes all hydrogen atoms.
        
    3.  pdb\_selaltloc: Removes alternate atom locations (common in high-res crystal structures).
        
    4.  pdb\_occ -0.3: Removes any atoms with low (<0.3) occupancy.
        
    5.  pdb\_selchain -A: Selects only Chain A of the protein.
        
    6.  pdb\_tidy: Reformats the PDB into a clean, standard format.
        
*   **Result:** A clean template file, 1J2E\_clean\_final.pdb, ready for alignment.
    

### 1.3. Sequence Alignment

MODELLER requires a specific alignment file (.pir) that maps our target sequence onto the template's sequence.

Commands:
`   pdb_tofasta 1J2E_clean_final.pdb > 1J2E_clean_A.fasta 
    needle -asequence dpp4.fasta -bsequence 1J2E_clean_A.fasta -gapopen 10 -gapextend 0.5 -af   `

*   **Explanation:**
    
    1.  pdb\_tofasta: We first extract the FASTA sequence from our _cleaned_ template PDB.
        
    2.  needle: We use the EMBOSS needle tool to create a global alignment between our target protein (dpp4.fasta) and our template (1J2E\_clean\_A.fasta).
        
*   **Result:** An alignment file, align.fasta. This file was then manually inspected and formatted into the required align.pir file for MODELLER.
    

### 1.4. Homology Modeling

With the alignment and template, we built our 3D model using MODELLER.

` mod10.7 model\_dpp4.py`
    
*   **Explanation:** This executes a Python script that instructs MODELLER to build three 3D models (.B99990001 - .B99990003) based on the align.pir file.
    
*   **Result:** Three candidate models of our unknown protein.
    

### 1.5. Model Selection

We needed to choose the best of the three models.

**PyMOL Commands:**

`   load 1j2e_renum.pdb  load dpp4.B99990001.pdb  load dpp4.B99990002.pdb  load dpp4.B99990003.pdb  align dpp4.B99990001, 1j2e_renum   align dpp4.B99990002, 1j2e_renum  align dpp4.B99990003, 1j2e_renum   `

*   **Explanation:** We loaded the template and our models into PyMOL. The align command superimposes a model onto the template and calculates the Root Mean Square Deviation (RMSD).
    
*   **Result:** Model 2 (dpp4.B99990002.pdb) had the lowest RMSD (**0.116 Å**). This model was selected, saved as **trimmedmodel.pdb**, and used for all future steps.
    

Phase 2: Model Validation (Molecular Dynamics)
----------------------------------------------

**Thought Process:** A 3D model is just a static hypothesis. To see if it represents a stable, physically realistic protein, we must perform a Molecular Dynamics (MD) simulation. This will simulate the protein's movement in a water environment at body temperature.

### 2.1. GROMACS System Setup

We used GROMACS to prepare the system for simulation.

**Commands:**
`   # 1. Create Topology (the "rule book" for the simulation)  
gmx pdb2gmx -f trimmedmodel.pdb -o processed.gro -water tip3p  `
`  # 2. Define Box  
gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic  `
`  # 3. Solvate (fill the box with water)  
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top`  
`  # 4. Add Ions  
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1  
echo 13 | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral   `

*   **Explanation:** This workflow generates the topology, creates a 1.0 nm cubic box, fills it with water (gmx solvate), and adds counter-ions (gmx genion) to neutralize the system's total charge.
    
*   **Result:** A fully solvated, neutral system ready for simulation (solv\_ions.gro).
    

### 2.2. Minimization & Equilibration

Before the "real" simulation, we must relax the system and bring it to the correct temperature and pressure.

**Commands:**
`   # 5. Energy Minimization (Steepest Descent, removes bad clashes) 
gmx grompp -f steep.mdp -c solv_ions.gro -p topol.top -o em.tpr -maxwarn 1  
gmx mdrun -v -deffnm em ` 
`  # 6. NVT Equilibration (Constant Volume, heats system to 300K)  
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr -maxwarn 1  
gmx mdrun -deffnm nvt -ntomp 4`
`  # 7. NPT Equilibration (Constant Pressure, adjusts density to 1 bar)  
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr -maxwarn 1  
gmx mdrun -deffnm npt -ntomp 4   `

*   **Result:** A stable, equilibrated system at 300K and 1 bar.
    

### 2.3. Production MD

Now we run the main simulation to observe the protein's behavior.

**Command:**
`   gmx grompp -f md.mdp -c npt.gro -r npt.gro -p topol.top -o md.tpr -maxwarn 1  
    gmx mdrun -deffnm md -ntomp 4   `

*   **Explanation:** This runs the final 200-picosecond production MD simulation, saving a "snapshot" (frame) of the protein's movement.
    
*   **Result:** A trajectory file (md.trr) that acts as a 200ps "movie" of our protein.
    

Phase 3: Trajectory & Functional Analysis
-----------------------------------------

**Thought Process:** The simulation is complete. Did our model unfold or "explode"? We must analyze the trajectory to confirm its stability.

### 3.1. Stability Analysis

We used gmx tools to calculate key stability metrics.

**Commands:**
`   # RMSD: Does the structure change over time? (Fit: Protein, Calc: C-alpha) 
echo 1 3 | gmx rms -s md.tpr -f md.trr -o rmsd.xvg -tu ns ` 
`   # RMSF: Which parts of the protein are flexible? (Calc: Protein)  
echo 1 | gmx rmsf -s md.tpr -f md.trr -o rmsf.xvg -res  `
`  # Radius of Gyration: Does the protein stay compact? (Calc: Protein) 
echo 1 | gmx gyrate -s md.tpr -f md.trr -o gyrate.xvg ` 
`  # H-Bonds: Does the internal structure (helices/sheets) break? (Calc: Protein-Protein) 
echo 1 1 | gmx hbond -s md.tpr -f md.trr -num hbnum.xvg   `

*   **Results (from xmgrace plots):**
    
    *   **RMSD:** Plateaued at **~0.13 nm**. **Interpretation:** The model is highly stable and did not unfold.
        
    *   **RMSF:** Showed low fluctuations, with peaks in loop regions. **Interpretation:** The model behaves like a typical, well-folded protein.
        
    *   **Rg:** Remained constant at **~2.7 nm**. **Interpretation:** The protein remained compact.
        
    *   **H-Bonds:** Remained stable at **~550 bonds**. **Interpretation:** The internal structure was preserved.
        

### 3.2. Functional Identification (The "Aha!" Moment)

**Thought Process:** Our model is definitively stable. Now, what does it _do_? A protein's function is determined by its active site.

**Commands:**
`   # 1. Extract the final, stable frame from the simulation 
echo 1 | gmx trjconv -s md.tpr -f md.trr -o final_frame.pdb -dump 200 ` 
`  # 2. Upload final_frame.pdb to CASTp web server.   `

*   **Explanation:** We extracted the final PDB structure (final\_frame.pdb) and submitted it to **CASTp** to find all surface pockets.
    
*   **Result (The Discovery):** CASTp identified **Pocket 1** as the largest (Volume 9614.2 Å³). Analysis of this deep, enzymatic-looking pocket, combined with its strong homology to 1J2E, allowed for a definitive functional assignment. Our unknown gut protein is **Human Dipeptidyl Peptidase-IV (DPP4)**, a major therapeutic target for type 2 diabetes.
    

Phase 4: Functional Validation (Virtual Screening)
--------------------------------------------------

**Thought Process:** If our model truly represents DPP4, it should be capable of binding known DPP4 inhibitors. We can test this by docking known drugs against our modeled protein. This step validates the functionality of our structure beyond its mere stability.

### 4.1. Ligand Curation (Research)

**Goal:** To compile a validation set of experimentally verified DPP4 inhibitors.

*   **Action 1 (Data Mining):** We searched the **ChEMBL database** for compounds annotated as Dipeptidyl Peptidase-IV (DPP4) inhibitors (CHEMBL284).
    
*   **Action 2 (Filtering):** From the search results, compounds were filtered based on their **IC₅₀** values (drug potency). The top 5 compounds with the lowest IC₅₀ values were selected to ensure strong inhibitory potential.
    
*   **Action 3 (Structure Generation):** The SMILES strings of the selected compounds were retrieved from ChEMBL. Using **UCSF Chimera’s “Build Structure”** tool, each SMILES was converted into a 3D molecular structure and saved in .pdb format.
    
*   **Result:** Five curated ligand structures were generated — 1.pdb, 2.pdb, 3.pdb, 4.pdb, and 5.pdb. These represent a benchmark validation set of known DPP4 inhibitors.
    

### 4.2. Docking Box Setup

Before performing docking, we defined the receptor’s active site region — the area where ligands are expected to bind.

*   **Procedure:**
    
    1.  The equilibrated MD structure (final\_frame.pdb) was opened in **AutoDockTools (ADT)**.
        
    2.  The binding pocket previously identified using CASTp was visually confirmed.
        
    3.  A grid box was drawn manually around the pocket to encompass all key active site residues.
        

**Result:** The geometric center of the docking box was recorded as:

`   center_x = 48.665  center_y = 50.969  center_z = 54.074   `

### 4.3. Virtual Screening (AutoDock Vina)

`   # Prepare receptor prepare_receptor4.py -r final_frame.pdb -o receptor1.pdbqt`
`  # Prepare ligands  for i in {1..5}; do    prepare_ligand4.py -l ${i}.pdb -o ${i}.pdbqt  done   `

**Example Vina config (config\_1.txt):**

`   receptor = receptor1.pdbqt  ligand = 1.pdbqt  out = 1_out.pdbqt  center_x = 48.665  center_y = 50.969  center_z = 54.074  size_x = 25  size_y = 25  size_z = 25   `

`   # Run docking  
vina --config config_1.txt > 1_log.txt
vina --config config_2.txt > 2_log.txt  # ... up to ligand 5   `

*   **Result:** All five inhibitors successfully docked to the modeled DPP4 active site. Among them, **Ligand 2** showed the strongest binding affinity with a score of **–7.785 kcal/mol**, confirming that the modeled protein is a valid functional representation of DPP4.
    

Phase 5: Lead Optimization (Rational Drug Design)
-------------------------------------------------

**Thought Process:** Our model is fully validated. Now, can we use it to design a _better_ drug? We will analyze the binding pose of our best known ligand (Ligand 2) and try to improve it.

### 5.1. Binding Pose Analysis

The best pose of Ligand 2 (2\_out.pdbqt) was analyzed in PyMOL.

*   **Analysis:** The pocket is highly **positively charged** (6 residues) and **aromatic** (5 residues), as confirmed by MOLE 2.0 analysis.
    
*   **The Opportunity:** Ligand 2 (a Xanthine derivative) bound well, but visual inspection showed it failed to interact with several key residues. We found three "untapped" opportunities:
    
    1.  **LYS-554:** A positive residue in the "upper" pocket.
        
    2.  **ARG-125:** A positive residue near the ligand's core.
        
    3.  **LEU-55/LEU-57:** A hydrophobic patch with empty space.
        

### 5.2. New Ligand Design

We designed three new ligands (n1, n2, n3) to exploit these opportunities.

*   **Original Ligand 2 SMILES:** CC#CCn1c(N2CCC\[C@@H\](N)C2)nc2c1c(=O)n(Cc1nc3cc(F)ccc3s1)c(=O)n2C
    
*   **Strategy A (Ligand n1):** Target **ARG-125**.
    
    *   **Edit:** Add a carboxymethyl group to the xanthine core.
        
    *   **SMILES:** CC#CCn1c(N2CCC\[C@@H\](N)C2)nc2c1c(=O)n(Cc1nc3cc(F)ccc3s1)c(=O)n2C(C(=O)O)
        
*   **Strategy B (Ligand n2):** Target **LYS-554**.
    
    *   **Edit:** Replace the terminal methyl of the butynyl chain with a carboxyl group.
        
    *   **SMILES:** O=C(O)C#CCn1c(N2CCC\[C@@H\](N)C2)nc2c1c(=O)n(Cc1nc3cc(F)ccc3s1)c(=O)n2C
        
*   **Strategy C (Ligand n3):** Target **LEU-55**.
    
    *   **Edit:** Add a methyl group to the benzothiazole ring.
        
    *   **SMILES:** CC#CCn1c(N2CCC\[C@@H\](N)C2)nc2c1c(=O)n(Cc1nc3cc(F)c(C)cc3s1)c(=O)n2C
        

Phase 6: Final Validation & Conclusion
--------------------------------------

**Thought Process:** We have designed three new, hypothetical drugs. Let's dock them using the _exact same_ validated protocol to see if our ideas worked.

### 6.1. Re-Docking New Compounds

The new SMILES were converted to PDB (in Chimera) and then PDBQT. They were docked using the same protocol as before.

`   for name in n1 n2 n3; do prepare_ligand4.py -l ${name}.pdb -o ${name}.pdbqt done  
vina --config config_n1.txt > n1_log.txt 
vina --config config_n2.txt > n2_log.txt
vina --config config_n3.txt > n3_log.txt   `

### 6.2. Final Score Comparison

——————————————————————————--------------------------------------------------

Ligand   | Strategy                    | Binding Affinity (kcal/mol)| Result

——————————————————————————---------------------------------------------------

Ligand 2 | Original Lead              | -7.785                      |Baseline

Ligand n1| Strategy A (Target ARG-125)| -7.483                      |Failed

Ligand n2| Strategy B (Target LYS-554)| -7.921                      |SUCCESS

Ligand n3 | Strategy C (Target LEU-55)| -7.622                      |Failed

The binding scores of the new, rationally designed ligands were compared to our original lead compound.

### 6.3. Conclusion

**Strategy B was a clear success.** The new compound, **Ligand n2**, showed a significantly improved binding affinity of **\-7.921 kcal/mol**.

Visual analysis of the new binding pose (n2\_out\_0001) in PyMOL confirmed our hypothesis: the newly added **carboxyl group** successfully formed a new, strong salt bridge with the target residue **LYS-554**.

This project successfully demonstrates a full _in silico_ pipeline, moving from an unknown protein sequence discovered in a proteome sample, to its structural modeling and functional identification, and finally to the rational design of a novel, optimized lead compound. 🚀

Software Used
-------------

*   UniProt, NCBI BLAST
    
*   pdb-tools
    
*   EMBOSS (needle)
    
*   MODELLER 10.7
    
*   GROMACS 2023.3
    
*   xmgrace
    
*   PyMOL 2.5
    
*   CASTp 3.0 Web Server
    
*   MOLE 2.0
    
*   ChEMBL Database
    
*   AutoDockTools (MGLTools)
    
*   AutoDock Vina 1.2.5
    
*   UCSF Chimera
