---
layout: post
title:  "MSY-based Reference Points"
date:   2014-03-23 06:50:11
categories: iscam overview
---

##MSY-based reference points

Fisheries references points based on the concept of maximum sustainable yield is fairly straight forwart to calculate for a single fishery.  However, when there are multiple fisheries simultaneously targeting the same stock, in addition to other non-target fisheries that may have a bycatch allocation, the calculus for determining optimal fishing mortality rate that maximizes yield is a bit more challenging.


In iSCAM, a great deal of effort has been spent on developing algorithms for estimating reference points for cases where there are two or more competeing gears that have different selectivity curves and fixed allocations associated with each gear. Much of the code for this work has been developed but more work is required on testing this code in more complex models.  The code is also written using C++ template classes and can be compiled directly in R using R to cpp.  The msy.cpp code is available at the following [Gist](https://gist.github.com/smartell/9695640)



{% highlight cpp %}
for(int i=1; i<=n; i++)
{
	f(i) = f(i) - df(i)/d2f(i);
}
{% endhighlight %}