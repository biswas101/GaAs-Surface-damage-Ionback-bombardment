# GaAs Surface Damage-Ionback-bombardment

This code can be used to find QE decay and power density on photocathode surface.<br/>
Ion back bombardment on the cathode surface can sputter away NEA material such as Cesium from acthode surface. Static Sputtering of Cesium can be calculated using SRIM Monte Carlo Simulation. <br/>
* Find QE map on the cathode surface<br/>
  * ```qe_map_nm.py``` can be used to find the ```Normalized QE``` and ```Ion Density```on the cathode surface.<br/>
  * ```qe_map_nm_2_f_2.py``` uses adaptive meshing to increase the performance of the calculation.<br/>
* FInd Power Density on Cathode Surface<br/>
  * ```pw_density_v2.py``` can be used to find the ```Power Density``` on Cathode Surface<br/>

## Prerequisites:

* Python 3
* or Pycharm with Proper Packages
* SRIM need to be installed for specific case of QE study. For Normalized Power density Estimation, SRIM is not required.

## Author

 * **Jyoti Biswas**

## License

This project is licensed under the ```MIT License``` - see the [LICENSE](LICENSE) file for details<br/>

## Acknowledgments

Great appreciation to **Dr. Edong Wang** for his valuable and constructive suggestions.

