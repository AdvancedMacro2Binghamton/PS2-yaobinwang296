# Problem Set 2

"Due" date (as in, be ready to discuss what you have done so far): Tuesday, Sep 19.

Yaobin's Comments:

My answers to question 1, 2, and 3 can be found in the m file named "VFI".

The plots for question 2 and 3 are saved in the folder named "plots".

Also in the last part of VFI.m, I simulated the economy with A_h = 1.1 and 
A_l = 0.678 and got a sequence of A, a sequence of K, and a sequence of Y.
Then I calculated the standard deviation of Y.

For calibrating A_h and A_l, I wrote my code in the m file called "Calibration". 
The basic logic of my calibration is based on the simulation I did in the VFI.m. 
I ran the calibration multiple times and the results imply that for A_h between 
1 to 1.0013 and for A_l between 0.9958 to 1, respectively, the standard deviation 
of output (Y) is around 2%, which is close to the 1.8% in the US data.