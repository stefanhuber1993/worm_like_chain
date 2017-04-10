# Worm Like Chain Simulator

## Description
Simulates a bundle of Worm-Like-Chains with a defined Persistence Length

## Dependencies
numpy, scipy, matplotlib

## Usage
```python
from WLC import ChainBundle
bundle = ChainBundle()
bundle.compute(n_chains=10, persistence_length=1, contour_length=10, n_segments=100)

bundle.data.shape
>>>(10, 100, 2) # 10 chains, containing 100 points, (x,y) coordinates

bundle.data
>>>array([[[ 0.        ,  0.        ],
        [ 0.1       ,  0.        ],
        [ 0.19048374, -0.04257573],
        ..., 
        [ 3.98511392,  5.23248311],
        [ 3.92634524,  5.31339197],
        [ 3.90761667,  5.41162252]],
       ...,        
       [[ 0.        ,  0.        ],
        [ 0.1       ,  0.        ],
        [ 0.19048374,  0.04257573],
        ..., 
        [-1.33431195,  1.46191013],
        [-1.30943588,  1.55876663],
        [-1.24568973,  1.63581485]]])

bundle.compute_persistence_length()
>>>1.0849898854583857

# Plotting
import matplotlib.pyplot as plt
bundle1 = ChainBundle()
bundle1.compute(n_chains=10, persistence_length=1, contour_length=10, n_segments=100)
bundle10 = ChainBundle()
bundle10.compute(n_chains=10, persistence_length=10, contour_length=10, n_segments=100)
bundle100 = ChainBundle()
bundle100.compute(n_chains=10, persistence_length=100, contour_length=10, n_segments=100)

bundle1.plot(color='red', label='Pers.Length 1 nm') 
bundle10.plot(color='green', label='Pers.Length 10 nm')
bundle100.plot(color='blue', label='Pers.Length 100 nm')

plt.legend()
plt.xlabel('x [nm]')
plt.ylabel('y [nm]')

plt.show()

```

## Source
Simulation procedure is explained in appendix of this paper:
* Castro, C. E., Su, H.-J., Marras, A. E., Zhou, L., & Johnson, J. (2015). Mechanical design of DNA nanostructures. Nanoscale, 7(14), 5913–5921. https://doi.org/10.1039/C4NR07153K_
