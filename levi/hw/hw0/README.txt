Hi Jessie! Sorry this is late! it wasn't too hard i just thought this was due 
next week + was occupied by AoT! Thanks for coming!!

I hope it follows the directions well - i wasn't sure if you wanted everything 
into one massive function because that felt a little bloated to me. Even doing
the linear fit and plotting in one function felt like pushing it a little -- a 
software dev friend of mine told me to have functions do ONLY one thing. if i was
wrong let me know!

To run it simply do:
python hw0.py --slopeval=3 --interceptval=0.5 --numpoints=4

NOTE: the numpoints argument is optional and defaults to 10, but i wanted to see
if the fit actually did worse with fewer points. After running it, you should get
some nice feedback with print statements and timings. See below for example:

EXAMPLE CMD OUTPUT:

points generated from input line: y = 4.0x+0.5
total time to generate points: 9.918212890625e-05 s
line fit to given points
time to fit line to data: 0.00013518333435058594 s
plotting beginning
time to plot: 0.15759015083312988 s
file [fitlinetopoints.png] saved in script's current directory
to finish the script, close the plot that has popped up.
plotting ending
total time for whole program: 2.572622060775757 s