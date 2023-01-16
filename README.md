<h1>Known Issues</h1>

* If the largest gap exists between the last prime in one thread, and first prime in second thread,
	then that gap will not be recognized.
* Efficiency issue: dividing N into equally sized chunks is not efficient, since it takes much longer for larger numbers to check if they are prime, than smaller numbers.
