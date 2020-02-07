# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 14:17:25 2020

@author: Libra
"""
import multiprocessing
import per_su_tca
import time

if __name__ == "__main__":
    process_list = []
    for epoc in ("LD", "ED", "DM"):
        for effect in ("impair", "improve"):
            for non_neg in (True, False):
                process = multiprocessing.Process(
                    target=per_su_tca.run_tca,
                    args=(160,),
                    kwargs={
                        "non_neg": non_neg,
                        "sep_blocks": False,
                        "epoc": epoc,
                        "effect": effect,
                    },
                )
                process.start()
                process_list.append(process)
                time.sleep(5)
    [process.join() for process in process_list]
