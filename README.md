MatRROMA is a MATLAB tool designed for the identification of dynamic properties of structures using the System Subspace Identification algorithm. It is cross-platform and runs on any operating system supported by MATLAB, including Windows, macOS, and GNU/Linux 64-bit.


This software performs modal identification using either single-setup (simultaneous) or multi-setup (fragmented) measurements. Regarding to multi-setup configuratins ,the software allows identification if there is at least one measurement point common to all sets (Reference-Based SSI) and also if there is at least one sensor common to two or more sets throughout the entire concatenation of sets (Roving-References SSI).  If you use MatRROMA please cite the following reference in your work (books, articles, reports, etc.):

Samuel González-Jiménez, Luis Alberto Padrón, Guillermo Manuel Álamo, Jacob Rodríguez Bordón. Operational Modal Analysis from multi-setup measurements with roving references and Stochastic Subspace Identification, _Engineering Structures_ 367, 123529, 2026. [Link](https://www.sciencedirect.com/science/article/pii/S0141029626014434).




## How to use?

1. **Command-Line Initialization:** The user runs the script via CLI, passing initial parameters such as the input data file, the number of block rows of the Hankel matrices, the system order, and other configuration settings.

2. **Interactive GUI & Pole Selection:** Once the initial data processing is complete, an interactive Graphical User Interface (GUI) is automatically launched. This interface allows the user to:
   - Define and adjust stabilization criteria.
   - Visualize the stabilization diagram.
   - Select the representative stable poles.

3. **Data Export:** Once the representative poles are selected, the tool asks the user to save those parameters along with other variables that help in the interpretation of the results.

**Additional Features:**
- Function to visualize the identified mode shapes.
- Merges data from multiple setups using the PoSER identification method.

For more details, please check the [User Manual](./docs/Manual.pdf). You can also find a complete identification example in the [Examples section](./docs/examples).

## License:

GPLv3, please see the [LICENSE](./LICENSE) file for details.

## Financing:

* **PID2023-151635OB-I00**, funded by the Ministerio de Ciencia, Innovación y Universidades and Agencia Estatal de Investigación (MICIU/AEI/10.13039/501100011033) of Spain, and FEDER/EU.

  <div align="center">
   <img src="./docs/img/MICIU+AEI.jpg" width="300">
  </div>
  
  
* **FPI2024010103**, research fellowship awarded to Samuel González-Jiménez from the Program of predoctoral fellowships from the Consejería de Universidades, Ciencia e Innovación y Cultura (Agencia Canaria de la Investigación, Innovación y Sociedad de la Información), Spain of the Gobierno de Canarias and Fondo Social Europeo Plus.

  <div align="center">
   <img src="./docs/img/Logo_ACIISI.jpg" width="300">
  </div>

