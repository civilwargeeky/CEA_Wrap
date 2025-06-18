import os


def getenv_t_f(env_variable:str, default:bool) -> bool:
  """
  Parse an environment variable for true or false values. Only looks at the first letter. Case insensitive. Values like "[y]es" or "[t]rue" or "1" are True
  and values like "[n]o" and "[f]alse" and "0" are False.

  :param env_variable: The environment variable to check
  :param default: The default value if the variable is not found or invalid
  :return: True or False
  """
  return {
    None: default,
    "t": True,
    "y": True,
    "1": True,
    "f": False,
    "n": False,
    "0": False,
  }.get(
    os.getenv(env_variable, "").lower()[:1],
    default
  )