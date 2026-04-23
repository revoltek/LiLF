from LiLF import lib_scheduler, lib_log, lib_walker


def run(step):
    log_dir = lib_log.Logger(f'pipeline-{step.kind}-{step.name}').log_dir
    logger = lib_log.logger
    s = lib_scheduler.Scheduler(dry_run=False, log_dir=log_dir)
    w = lib_walker.Walker(f'pipeline-{step.kind}-{step.name}.walker')

    logger.info(f'Started step: facetselfcal - {step.name}')
    w.alldone()
