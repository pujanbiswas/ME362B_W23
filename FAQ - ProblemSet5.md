# Bugs and Errors: Problem Set 5

## 1. In _Part 2_, which definition of IDT should I use?

The answer shouldn't matter much on the definition selected. I would recommend using 'T' just to be consistent with what is asked for in later questions.

## 2. My ignition delay times/sensitivity results do not agree at all between mechanisms? Is this right?

No, they should be agreeing pretty closely (within ~50%). If they are not, and especially if you are getting IDT values on the order of 1e-11, make sure you have implemented the bug fixes called out in the announcement on Canvas posted Feb. 24. 

## 3. In _Part 3_, my "most-negative reactions" don't agree between mechanisms. Am I doing something wrong?

No, I misread my own sensitivity plots when writing the problem. Sorry!

You may either ignore the part stating the reactions should be the same and characterize all the reactions that appear as the most-negatively sensitive for any mechanism, or you may consider the reaction "CH3 + HO2 <-> CH3O + OH" as the most negatively sensitive (which should be the case for 2 of the 3 mechanisms).

## 4. In _Part 4_, I get a Cantera error: Unknown Species 'C3H8'

There was an error in the problem statement in earlier versions of the Problem Set template; please use __methane__ as the fuel in part 4, as in the prior parts of the assignment.

The error is raised because the reaction mechanisms do not contain propane as a species.

*********

# Frequently Asked Questions: Problem Set 5

## 1. How do I save my work in _Binder_?

a. As a Jupyter notebook to continue working on later:
Use _File > Download as > Notebook (.ipynb)_ to save a copy of your modified notebook in its current state.

b. As a PDF for submission: 
Use _File > Print Preview_, then _Print_ the preview as a .PDF file.

## 2. How do I load my work to continue working?

Open a _Binder_ session and navigate to the repository view:

* Option 1: Click here: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ajsusa/me362b_winter2021/602b07cd40f984396f5634f443f8dfba903447a6), then open the _notebooks_ folder
* Option 2: Use the _Launch Binder_ link from the Problem Set assignment to open the problem set 5 template, then click _File > Open_

In the upper-right corner, click _Upload_, then navigate to and select your file.

Rename your file for upload if needed, then click the _Upload_ button in line with your file.

Click the file to open it in the interactive session.

**_Note:_** you have to upload your file into the _notebooks_ directory, else the filepaths to the mechanism files will break 

## 3. Can I run my notebook locally on my local computer, instead of through _Binder_?

Yes you can, but for the sake of this assignment, it's not recommended or worth the effort.

To run your notebook on your local computer, you would need:

* Python 3 installed
* A Python interpretter configured with Cantera, numpy, matplotlib, jupyter, and ipykernel installed
* Your interpretter registered as Jupyter kernel
* Resetting the notebook kernel for your assignment to your local interpretter

None of this is particularly difficult, but it introduced a good bit of extra time, a number of extra steps,
and the need to install software that can all be avoided by using the preconfigured _Binder_ environment. 
