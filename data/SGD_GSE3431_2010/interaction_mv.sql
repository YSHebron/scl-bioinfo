create materialized view interaction_mv
pctused 60 pctfree 10 tablespace VIEW2
storage ( pctincrease 0 initial 10M NEXT 10M ) 
build immediate
never refresh
enable query rewrite 
as
select  i.interaction_no, 
fa.feature_name feature_a, fb.feature_name feature_b,
fa.gene_name gene_a, fb.gene_name gene_b,
regexp_substr(an.headline, '[[:alnum:]]+[[:alnum:]|[:print:]|[:punct:]][^,|^;|^.]+', 1) headline,
i.experiment_type, i.interaction_type, b.action, i.source, 
replace (u.url, '_SUBSTITUTE_THIS_', d.dbxref_id) url, 
i.annotation_type,
(CASE when i.modification = 'No Modification'
      then ir.note
      else i.modification || '; ' || ir.note
      end ) notes,
(CASE when ph.qualifier is not null 
      then ph.observable || ':' || ph.qualifier
      else ph.observable
      end ) phenotype,
r.dbxref_id, r.pubmed, r.citation
from bud.interaction i, bud.feat_interact a, bud.feat_interact b,
bud.interact_pheno ip, bud.phenotype ph,
bud.reference r, bud.interact_ref ir, 
bud.feature fa, bud.feature fb, bud.feat_annotation an,
bud.dbxref_feat df, bud.dbxref d, bud.dbxref_url du, 
bud.url u, bud.web_display wd
where a.action != b.action
and fa.feature_no != fb.feature_no
and a.feature_no = fa.feature_no
and b.feature_no = fb.feature_no
and an.feature_no = fb.feature_no 
and ip.interaction_no(+) = i.interaction_no
and ip.phenotype_no      = ph.phenotype_no(+)
and i.interaction_no = a.interaction_no
and i.interaction_no = b.interaction_no
and r.reference_no = ir.reference_no
and ir.interaction_no = i.interaction_no
and df.feature_no = fb.feature_no
and df.dbxref_no = d.dbxref_no
and d.dbxref_no = du.dbxref_no
and du.url_no = u.url_no
and u.url_no = wd.url_no
and wd.label_name = i.source
and wd.web_page_name = 'Interaction'
union
select i.interaction_no, 
fa.feature_name feature_a, fb.feature_name feature_b,
fa.gene_name gene_a, fb.gene_name gene_b,
regexp_substr(an.headline, '[[:alnum:]]+[[:alnum:]|[:print:]|[:punct:]][^,|^;|^.]+', 1) headline,
i.experiment_type, i.interaction_type, 'Self', i.source, 
replace (u.url, '_SUBSTITUTE_THIS_', d.dbxref_id) url,
i.annotation_type, 
(CASE when i.modification = 'No Modification'
      then ir.note
      else i.modification || '; ' || ir.note
      end ) notes,
(CASE when ph.qualifier is not null 
      then ph.observable || ':' || ph.qualifier
      else ph.observable
      end ) phenotype,
r.dbxref_id, r.pubmed, r.citation
from bud.interaction i, bud.feat_interact a, bud.feat_interact b,
bud.interact_pheno ip, bud.phenotype ph,
bud.reference r, bud.interact_ref ir, 
bud.feature fa, bud.feature fb, bud.feat_annotation an,
bud.dbxref_feat df, bud.dbxref d, bud.dbxref_url du, 
bud.url u, bud.web_display wd
where a.action = 'Bait'
and  b.action = 'Hit'
and fa.feature_no = fb.feature_no
and a.feature_no = fa.feature_no
and b.feature_no = fb.feature_no
and an.feature_no = fb.feature_no
and ip.interaction_no(+) = i.interaction_no
and ip.phenotype_no      = ph.phenotype_no(+)
and i.interaction_no = a.interaction_no
and i.interaction_no = b.interaction_no
and r.reference_no = ir.reference_no
and ir.interaction_no = i.interaction_no
and df.feature_no = fb.feature_no
and df.dbxref_no = d.dbxref_no
and d.dbxref_no = du.dbxref_no
and du.url_no = u.url_no
and u.url_no = wd.url_no
and wd.label_name = i.source
and wd.web_page_name = 'Interaction';

grant select on interaction_mv to curator
/

grant select on interaction_mv to dbselect
/

