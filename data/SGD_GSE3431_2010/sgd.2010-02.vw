-- C:\Documents and Settings\Gail Binkley\Desktop\tables\bud2.2010-02.vw
--
-- Generated for Oracle 10g on Fri Feb 12  09:33:21 2010 by Server Generator 10.1.2.9.28
 

PROMPT Creating View 'SAGE_CONDITIONS'
CREATE OR REPLACE FORCE VIEW SAGE_CONDITIONS
 AS SELECT se.sage_tag_seq, sum(decode(condition, 'L', exp_value, NULL)) L,
sum (decode(condition, 'S', exp_value, NULL)) S,
sum (decode(condition, 'G2M', exp_value, NULL)) G2M
from sage_expression se
group by se.sage_tag_seq
/

