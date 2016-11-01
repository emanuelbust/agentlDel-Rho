PYTHONS = $(wildcard *.py)
HELPS = $(wildcard *.help)
BASHS = $(wildcard *.sh)
BACKUP = simBackup.tar.gz

$(BACKUP): $(PYTHONS) $(HELPS) $(BASHS)
	tar -cvzf $(BACKUP) *.sh *.py *.help Makefile

.PHONY: clean
clean:
	rm -f $(BACKUP)
