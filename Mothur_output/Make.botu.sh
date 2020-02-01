sed s/Otu/Botu/g 16s_bark_taxonomy_fa19_AJO.cons.taxonomy > 16s_bark_taxonomy_fa19_AJO.cons.taxonomy.botu
sed s/Otu/Botu/g 16s_soil_taxonomy_fa19_AJO.cons.taxonomy > 16s_soil_taxonomy_fa19_AJO.cons.taxonomy.botu
sed s/Otu/Botu/g 16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared > 16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared.botu
sed s/Otu/Botu/g 16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared > 16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared.botu

sed s/\"//g 16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared.botu | sed $'s/BAC\t/BAC\t1775\t/g' | sed $'1!s/^/0.03\t/g'  > temp && mv temp 16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared.botu
sed $'1s/^/label\tGroup\tnumBotus\t/' 16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared.botu > temp && mv temp 16s_bark_otutable_rarefied_nosingleton_fa19_ajo.shared.botu
sed s/\"//g 16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared.botu | sed $'s/BAC\t/BAC\t21433\t/g'| sed $'1!s/^/0.03\t/g'  > temp && mv temp 16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared.botu
sed $'1s/^/label\tGroup\tnumBotus\t/' 16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared.botu > temp && mv temp 16s_soil_otutable_rarefied_nosingleton_fa19_ajo.shared.botu
