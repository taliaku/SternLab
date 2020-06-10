#! python/python-anaconda3.2019.7

def check_queue(queue):
    allowed_queues = ["inf", "hugemem", "pup-interactive", "parallel", "adis", "adis-long", "tzachi@power9"] 
    if queue not in allowed_queues:
        raise Exception(f"Sorry but queue must be one of {allowed_queues}, not '{queue}'")
