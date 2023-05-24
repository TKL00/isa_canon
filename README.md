# Individuel Studieaktivitet - Graph Canonization
af Kasper Halkjær Beider og Tobias Klink Lehn, datalogistuderende ved Institut for Matematik og Datalogi på Syddansk Universitet, Odense.


Bibliotektets implementering af nauty's grafkanoniseringsalgoritme ligger i src/canonization.py. Funktionen _graph_canon_ beregner den kanoniske repræsentation af en inputgraf. Der er mulighed for at implementere sin egen "target cell selector" og give denne som parameter til funktionen. Desuden kan man vælge om "traces" skal slås til - dette kan påvirke køretiden positivt for nogle inputgrafer.

