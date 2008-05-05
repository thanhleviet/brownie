#include "gport.h"

#include <iostream>

int main ()
{
	GPostscriptPort g;

    g.StartPicture ("test.ps");
    g.DrawLine (10, 10, 100, 20);

    GBaseFont f;

    g.SetCurrentFont (f);
    g.DrawText (20, 20, "Hello");

    GRect r (0, 300, 250, 320);

    g.DrawRect (r);

    g.EndPicture ();
    
	cout << endl << endl << "Done...press any key to end";

	char ch;
	cin >> ch;

	return 0;



}